"""OU2024 recovery in bins of true-host r magnitude, and the OU2024 / DP0.2 comparison at fixed host
magnitude.

OU2024 hosts are drawn ~1 dex too low in stellar mass (OU2024 paper 3.2.4), so they are fainter than
DC2's and the population-averaged recovery is not comparable between tiers; at fixed host magnitude
the association is the same geometry and photometry problem in both.

  by_host_mag     json: SBW / DLR recovery, SBW and DLR hostless fractions and SBW-DLR agreement per
                  host r bin, for all transients and for z < 1, plus the r < 24 summaries.
  dp02_comparison figure: SBW and DLR recovery versus host r for both tiers (OU2024 all z and z < 1)
                  and the host-magnitude distributions, with OU2024 reweighted to DP0.2's host r
                  distribution.
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dp02_fig import binned_fraction, load, subtitle
from dp02_lib import COLORS, HOSTLESS_CUT, sbw_hostless, wilson

MAG_EDGES = [10, 20, 22, 23, 24, 25, 26, 27, 35]
FINE = np.arange(16, 29.01, 0.5)


def stats(d, cut):
    n = len(d)
    hl = sbw_hostless(d, cut)
    ps, ls, hs = wilson(int(d.sbw_correct.sum()), n)
    pd_, ld, hd = wilson(int(d.dlr_correct.sum()), n)
    return {"n": n, "sbw_recovery": ps, "sbw_recovery_ci68": [ls, hs], "dlr_recovery": pd_, "dlr_recovery_ci68": [ld, hd],
            "sbw_recovery_with_hostless_cut": float((d.sbw_correct.to_numpy() & ~hl).mean()) if n else np.nan,
            "sbw_hostless_fraction": float(hl.mean()) if n else np.nan,
            "dlr_hostless_fraction": float(d.dlr_hostless.mean()) if n else np.nan,
            "agreement": float(d.agree.mean()) if n else np.nan,
            "true_host_in_candidates": float(d.true_in_cand.mean()) if n else np.nan,
            "dlr_right_sbw_wrong": int((d.dlr_correct & ~d.sbw_correct).sum()),
            "sbw_right_dlr_wrong_or_none": int((d.sbw_correct & ~d.dlr_correct).sum())}


def by_host_mag(df, cut, out):
    import pandas as pd
    res = {"bins_r": MAG_EDGES, "sbw_hostless_cut": cut}
    for name, d in [("all_z", df), ("z_lt_1", df[df.sn_redshift < 1])]:
        b = pd.cut(d.host_mag_r, MAG_EDGES)
        res[name] = {"all": stats(d, cut), "host_r_lt_24": stats(d[d.host_mag_r < 24], cut),
                     "bins": {str(k): stats(g, cut) for k, g in d.groupby(b, observed=False)}}
    res["settings"] = df.attrs
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def reweighted(ou, dp, flag):
    """OU2024 recovery reweighted to DP0.2's host-r distribution (bins of FINE; empty OU2024 bins dropped)."""
    k_ou, _ = np.histogram(ou.host_mag_r[ou[flag]], FINE)
    n_ou, _ = np.histogram(ou.host_mag_r, FINE)
    n_dp, _ = np.histogram(dp.host_mag_r, FINE)
    ok = n_ou > 0
    w = n_dp[ok] / n_dp[ok].sum()
    return float((w * k_ou[ok] / n_ou[ok]).sum()), float(n_dp[~ok].sum() / n_dp.sum())


def dp02_comparison(ou, dp, out):
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), width_ratios=[1, 1, 0.9])
    for ax, (oname, o) in zip(axes[:2], [("OU2024, all z", ou), ("OU2024, z < 1", ou[ou.sn_redshift < 1])]):
        for d, tier, ls in [(dp, "DP0.2", "-"), (o, oname, "--")]:
            for flag, meth, color in [("sbw_correct", "SBW", COLORS["sbw"]), ("dlr_correct", "DLR", COLORS["dlr"])]:
                c, p, lo, hi, n = binned_fraction(d.host_mag_r.to_numpy(float), d[flag].to_numpy(), FINE)
                ok = n >= 30
                ax.fill_between(c[ok], lo[ok], hi[ok], color=color, alpha=0.12, lw=0)
                ax.plot(c[ok], p[ok], color=color, ls=ls, label=f"{meth}, {tier}")
        rs, miss = reweighted(o, dp, "sbw_correct")
        rd, _ = reweighted(o, dp, "dlr_correct")
        b24 = lambda d, f: 100 * d[f][d.host_mag_r < 24].mean()
        ax.text(0.02, 0.03, f"{oname} reweighted to DP0.2 host r: SBW {100 * rs:.1f}%, DLR {100 * rd:.1f}%\n"
                            f"host r < 24: DP0.2 SBW {b24(dp, 'sbw_correct'):.1f} / DLR {b24(dp, 'dlr_correct'):.1f}%, "
                            f"OU2024 {b24(o, 'sbw_correct'):.1f} / {b24(o, 'dlr_correct'):.1f}%",
                transform=ax.transAxes, fontsize=7, va="bottom")
        ax.set_ylim(0, 1.02); ax.set_xlabel("true host $r$ [mag]"); ax.set_ylabel("fraction top host = true host")
        ax.set_title(f"DP0.2 (N = {len(dp):,}) vs {oname} (N = {len(o):,})", fontsize=8)
        ax.legend(fontsize=7, loc="center left")
    ah = axes[2]
    for d, lab, color in [(dp, "DP0.2", "#333333"), (ou, "OU2024, all z", "#7b4ea3"),
                          (ou[ou.sn_redshift < 1], "OU2024, z < 1", "#b07cc6")]:
        ah.hist(d.host_mag_r, bins=FINE, density=True, histtype="step", color=color, lw=1.5,
                label=f"{lab} (median {np.nanmedian(d.host_mag_r):.1f})")
    ah.set_xlabel("true host $r$ [mag]"); ah.set_ylabel("density"); ah.legend(fontsize=7, loc="upper left")
    ah.set_title("host magnitude distributions (sample)", fontsize=8)
    fig.suptitle("Host recovery at fixed host magnitude, DP0.2 vs OU2024\nOU2024: " + subtitle(ou)
                 + "  |  DP0.2: " + subtitle(dp), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["by_host_mag", "dp02_comparison"])
    ap.add_argument("--associations", required=True, help="ou24_associations")
    ap.add_argument("--dp02-associations")
    ap.add_argument("--hostless-cut", choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = load(a.associations)
    if a.kind == "by_host_mag":
        by_host_mag(df, HOSTLESS_CUT[a.hostless_cut], a.output)
    else:
        dp02_comparison(df, load(a.dp02_associations), a.output)


if __name__ == "__main__":
    main()
