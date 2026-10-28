"""Do OU2024 transients sit on the rendered galaxy or on its mirror image (the DC2 bug)?

For elongated hosts (skyCatalog q < 0.5 of the dominant component, B/T_i split) and transients beyond
1 R_e of that component (elliptical offset, ou24_sample.host_offsets), the SN offset angle psi (from +x
West toward +y North around the host truth position) is compared with
* the skyCatalog major-axis angle beta, as psi - beta (direct) and psi + beta (mirrored), and
* the PSF-subtracted r-band moment angle phi of the matched object (nearest within 1 arcsec, both
  matchers agreeing, moment ellipticity > 0.1), as psi - phi and psi + phi.
The moment test is model-free: phi is what the pipeline measured on the rendered image, so
<cos 2(psi - phi)> > 0 with <cos 2(psi + phi)> ~ 0 means the transients follow the rendered galaxy,
and the reverse means they follow its mirror image (DC2: 0.005 against 0.59 for q < 0.3).
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from dp02_fig import plt as _plt  # noqa: F401  (rcParams)
from dp02_lib import BULGE_DOMINATED_BT, COLORS, moment_shape
from ou24_sample import host_offsets
from sbw_host.core import sky_to_plane

INT64 = {pa.int64(): pd.Int64Dtype()}.get
Q_MAX, R_MIN, E_MIN = 0.5, 1.0, 0.1


def wrap(d):
    """Axis angle difference (radians) folded to [-90, 90) degrees."""
    return np.degrees((d + np.pi / 2) % np.pi - np.pi / 2)


def cos2(d):
    c = np.cos(2 * np.deg2rad(d))
    return {"mean": float(c.mean()) if len(c) else np.nan,
            "err": float(c.std(ddof=1) / np.sqrt(len(c))) if len(c) > 1 else np.nan, "n": int(len(c))}


def load(sn_table, match):
    cols = ["sn_id", "sn_ra", "sn_dec", "class", "host_id", "host_ra", "host_dec", "host_mag_r", "host_hlr_bulge",
            "host_hlr_disk", "host_q_bulge", "host_q_disk", "host_beta_bulge", "host_beta_disk", "host_bt_i", "a_rot"]
    d = pq.read_table(sn_table, columns=cols).to_pandas(ignore_metadata=True, types_mapper=INT64)
    m = pq.read_table(match).to_pandas(ignore_metadata=True, types_mapper=INT64)
    m = m[[c for c in m.columns if c.startswith("nearest_") or c in ("host_id", "matchers_agree", "mag_objectId")]]
    d = d[(d["class"] != "FIXMAG") & d.host_ra.notna()].merge(m, on="host_id", how="left")
    off = host_offsets(d)
    bd = off["bulge_dominated"]
    d["bulge_dominated"] = bd
    d["offset_re"] = off["own"]
    d["q"] = np.where(bd, d.host_q_bulge, d.host_q_disk).astype(float)
    d["beta"] = np.deg2rad(np.where(bd, d.host_beta_bulge, d.host_beta_disk).astype(float))
    x, y = sky_to_plane(d.sn_ra.to_numpy(float), d.sn_dec.to_numpy(float), d.host_ra.to_numpy(float),
                        d.host_dec.to_numpy(float))
    d["psi"] = np.arctan2(y, x)
    g = lambda k: d[f"nearest_r_{k}"].to_numpy(float)
    q, phi, l1, l2, valid = moment_shape(g("ixx"), g("iyy"), g("ixy"), g("ixxPSF"), g("iyyPSF"), g("ixyPSF"))
    with np.errstate(invalid="ignore", divide="ignore"):
        e = (l1 - l2) / (l1 + l2)
    d["phi"] = np.deg2rad(phi)
    d["mom_ok"] = valid & (e > E_MIN) & d.nearest_objectId.notna().to_numpy() & d.matchers_agree.fillna(False).to_numpy(bool)
    d["sel"] = (d.q < Q_MAX) & (d.offset_re > R_MIN) & np.isfinite(d.beta)
    return d


def metrics(d):
    s = d[d.sel]
    sm = s[s.mom_ok]
    out = {"selection": f"dominant-component q < {Q_MAX} (B/T_i > {BULGE_DOMINATED_BT} -> bulge), elliptical offset > "
                        f"{R_MIN} R_e; moments: nearest object within 1 arcsec, both matchers agree, e > {E_MIN}",
           "psi_minus_beta": cos2(wrap(s.psi - s.beta)), "psi_plus_beta": cos2(wrap(s.psi + s.beta)),
           "psi_minus_arot": cos2(wrap(s.psi - np.deg2rad(s.a_rot))),
           "psi_minus_phi": cos2(wrap(sm.psi - sm.phi)), "psi_plus_phi": cos2(wrap(sm.psi + sm.phi)),
           "beta_minus_phi": cos2(wrap(sm.beta - sm.phi)), "beta_plus_phi": cos2(wrap(sm.beta + sm.phi))}
    out["split"] = {}
    for name, m in [("q_lt_0p3", s.q < 0.3), ("bulge_dominated", s.bulge_dominated), ("disk_dominated", ~s.bulge_dominated)]:
        ss = s[m]; ssm = ss[ss.mom_ok]
        out["split"][name] = {"psi_minus_beta": cos2(wrap(ss.psi - ss.beta)), "psi_plus_beta": cos2(wrap(ss.psi + ss.beta)),
                              "psi_minus_phi": cos2(wrap(ssm.psi - ssm.phi)), "psi_plus_phi": cos2(wrap(ssm.psi + ssm.phi))}
    direct, mirror = out["psi_minus_phi"], out["psi_plus_phi"]
    sig = lambda r: r["mean"] / r["err"]
    if sig(direct) > 5 and sig(direct) > 3 * abs(sig(mirror)):
        v = "rendered: transients follow the galaxy as rendered (no DC2-like mirror bug)"
    elif sig(mirror) > 5 and sig(mirror) > 3 * abs(sig(direct)):
        v = "mirror: transients follow the mirror image of the rendered galaxy, as in DC2"
    else:
        v = "inconclusive"
    out["verdict"] = v
    out["dc2_reference"] = "DC2 (q < 0.3): <cos 2(psi - PA_true)> = 0.59, against the rendered (moment) axis 0.005"
    return out


def figure(d, res, out):
    s = d[d.sel]; sm = s[s.mom_ok]
    edges = np.linspace(-90, 90, 37)
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.6))
    for ax, (a_, b_, la, lb, n_, title) in zip(axes.flat[:3], [
            (wrap(s.psi - s.beta), wrap(s.psi + s.beta), "$\\psi - \\beta$", "$\\psi + \\beta$ (mirror)", len(s),
             "(a) SN angle vs skyCatalog major axis $\\beta$"),
            (wrap(sm.psi - sm.phi), wrap(sm.psi + sm.phi), "$\\psi - \\phi$", "$\\psi + \\phi$ (mirror)", len(sm),
             "(b) SN angle vs pipeline moment axis $\\phi$"),
            (wrap(sm.beta - sm.phi), wrap(sm.beta + sm.phi), "$\\beta - \\phi$", "$\\beta + \\phi$", len(sm),
             "(c) skyCatalog vs pipeline axis")]):
        ax.hist(a_, bins=edges, density=True, histtype="step", color=COLORS["sbw"], lw=1.5,
                label=f"{la}: $\\langle\\cos 2\\Delta\\rangle$ = {np.cos(2 * np.deg2rad(a_)).mean():.3f}")
        ax.hist(b_, bins=edges, density=True, histtype="step", color=COLORS["dlr"], lw=1.5,
                label=f"{lb}: $\\langle\\cos 2\\Delta\\rangle$ = {np.cos(2 * np.deg2rad(b_)).mean():.3f}")
        ax.axhline(1 / 180, color=COLORS["ref"], ls=":", lw=1)
        ax.set_xlabel("angle difference [deg]"); ax.set_ylabel("density"); ax.set_xlim(-90, 90)
        ax.set_title(f"{title}, N = {n_:,}", fontsize=8, loc="left"); ax.legend(fontsize=7, loc="upper right")
    ax = axes[1, 1]
    qe = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
    base = d[(d.offset_re > R_MIN) & np.isfinite(d.beta)]
    qc = 0.5 * (qe[1:] + qe[:-1])
    for col, lab, color, ls, mom in [("psi_minus_beta", "$\\psi - \\beta$", COLORS["sbw"], "-", False),
                                     ("psi_plus_beta", "$\\psi + \\beta$", COLORS["dlr"], "-", False),
                                     ("psi_minus_phi", "$\\psi - \\phi$", COLORS["sbw"], "--", True),
                                     ("psi_plus_phi", "$\\psi + \\phi$", COLORS["dlr"], "--", True)]:
        y, e = [], []
        for lo, hi in zip(qe[:-1], qe[1:]):
            b = base[(base.q >= lo) & (base.q < hi)]
            if mom:
                b = b[b.mom_ok]
            ang = {"psi_minus_beta": b.psi - b.beta, "psi_plus_beta": b.psi + b.beta,
                   "psi_minus_phi": b.psi - b.phi, "psi_plus_phi": b.psi + b.phi}[col]
            r = cos2(wrap(ang)); y.append(r["mean"]); e.append(r["err"])
        ax.errorbar(qc, y, yerr=e, color=color, ls=ls, marker="o", ms=3, capsize=2, lw=1, label=lab)
    ax.axhline(0, color=COLORS["ref"], lw=1)
    ax.set_xlabel("skyCatalog axis ratio q (dominant component)"); ax.set_ylabel("$\\langle\\cos 2\\Delta\\rangle$")
    ax.set_title(f"(d) alignment vs q, offset > {R_MIN:g} $R_e$", fontsize=8, loc="left"); ax.legend(fontsize=7)
    fig.suptitle(f"OU2024 SN placement vs host orientation: {res['verdict']}\n"
                 f"q < {Q_MAX}, offset > {R_MIN:g} $R_e$; $\\psi$, $\\beta$, $\\phi$ from +x (West) toward +y (North)",
                 fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["figure", "metrics"])
    ap.add_argument("--sn-table", required=True)
    ap.add_argument("--match", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    d = load(a.sn_table, a.match)
    res = metrics(d)
    if a.kind == "metrics":
        with open(a.output, "w") as f:
            json.dump(res, f, indent=2, default=float)
    else:
        figure(d, res, a.output)


if __name__ == "__main__":
    main()
