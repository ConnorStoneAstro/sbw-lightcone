"""SN offset from the truth host in units of the host's truth half-light radius (bulge, disk, and the
dominant component), against the SBW / DLR failure fractions, split into bulge- and disk-dominated hosts."""
from __future__ import annotations

import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq

from dp02_fig import binned_fraction, load, subtitle
from dp02_lib import BULGE_DOMINATED_BT, COLORS, HOSTLESS_CUT, host_offsets, sbw_hostless

EDGES = np.geomspace(0.03, 100, 29)
XLAB = {"bulge": "offset / bulge $R_e$", "disk": "offset / disk $R_e$",
        "own": "offset / $R_e$ of dominant component"}


def failure_flags(df, cut):
    s, d, t = df.sbw_objectId.to_numpy(), df.dlr_objectId.to_numpy(), df.host_objectId.to_numpy()
    return {"SBW top host wrong": (COLORS["sbw"], "-", s != t),
            "DLR host wrong or none": (COLORS["dlr"], "-", d != t),
            "both wrong, same object": ("#222222", "-", (s == d) & (s >= 0) & (s != t)),
            f"SBW hostless (flb < {cut:g})": (COLORS["sbw"], ":", sbw_hostless(df, cut))}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--shapes", required=True)
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    cut = HOSTLESS_CUT[a.hostless_cut]
    df = load(a.associations)
    attrs = df.attrs
    sh = pq.read_table(a.shapes).to_pandas(ignore_metadata=True).rename(columns={"galaxy_id": "host_galaxy"})
    df = df.merge(sh.drop(columns=["ra", "dec", "redshift"]), on="host_galaxy", how="inner")
    df.attrs = attrs
    off = host_offsets(df)
    flags = failure_flags(df, cut)
    bd = off["bulge_dominated"]
    classes = [(f"B/T$_i$ > {BULGE_DOMINATED_BT:g}", bd), (f"B/T$_i$ $\\leq$ {BULGE_DOMINATED_BT:g}", ~bd)]

    fig, axes = plt.subplots(3, 3, figsize=(10.5, 8.4), sharex=True, height_ratios=[3, 3, 2])
    summary = {}
    for j, comp in enumerate(("bulge", "disk", "own")):
        x = off[comp]
        for i, (cname, sel) in enumerate(classes):
            ax = axes[i, j]
            for lab, (color, ls, f) in flags.items():
                c, p, lo, hi, n = binned_fraction(x[sel], f[sel], EDGES)
                ok = n >= 20
                ax.fill_between(c[ok], lo[ok], hi[ok], color=color, alpha=0.12, lw=0)
                ax.plot(c[ok], p[ok], color=color, ls=ls, label=lab)
                if comp == "own":
                    summary.setdefault(cname, {})[lab] = {"x": c[ok].round(3).tolist(), "frac": p[ok].round(4).tolist(),
                                                          "n": n[ok].tolist()}
            ax.set_xscale("log"); ax.set_ylim(0, 1.02)
            ax.set_title(f"{cname} (N = {int(sel.sum()):,})", fontsize=8)
            if j == 0:
                ax.set_ylabel("fraction of SNe")
        ah = axes[2, j]
        for (cname, sel), color in zip(classes, ["#8c2d04", "#2c7fb8"]):
            ah.hist(x[sel & np.isfinite(x)], bins=EDGES, histtype="step", color=color, label=cname)
        ah.set_yscale("log"); ah.set_xlabel(XLAB[comp])
        if j == 0:
            ah.set_ylabel("SNe")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    axes[2, 0].legend(fontsize=7, loc="upper right")
    fig.suptitle("Failures versus SN offset in host truth half-light radii (elliptical, cosmoDC2 shapes)\n"
                 + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(a.output); plt.close(fig)
    json.dump(summary, sys.stderr, indent=None)
    print(file=sys.stderr)


if __name__ == "__main__":
    main()
