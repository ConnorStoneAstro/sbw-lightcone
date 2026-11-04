"""Run the shared DP0.2 / DP2 scripts on the DES-SN5YR tier, and the DES-only recovery metric.

    python src/des_run.py associate <src/dp02_associate.py arguments>
    python src/des_run.py sweep     <src/dp02_sweep.py arguments>     (bands griz)
    python src/des_run.py fig       <src/dp02_fig.py arguments>       (--tier des)
    python src/des_run.py dp2fig    <src/dp2_fig.py arguments>        (bands griz)
    python src/des_run.py recovery  --associations ... --random ... --sample ... --sweep ... --hostless-cut ... --p-threshold ... --output ...

``des_lib.install()`` swaps the DES light models into ``dp02_lib``; the band loops of the sweep and DP2 figures
are set to griz (DES DR2 has no u band, and Y is shallow and unused by DES-SN). On DES the reference is the
published DES-SN host, a DLR assignment on the deep SN-field coadds, not truth: every figure made through
``fig`` / ``dp2fig`` has its wording rewritten at save time ("true host" -> "published host", "DP2" ->
"DES-SN5YR"), so the shared figure code stays untouched.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import des_lib

des_lib.install()

import dp02_lib as L   # noqa: E402
from dp02_lib import HOSTLESS_CUT, P_TRUE, dlr_cols, prob_rows, sbw_hostless, sbw_success, wilson   # noqa: E402

SWEEP_BANDS = "griz"
WORDS = [("true host", "published host"), ("= true", "= published"),
         ("DLR host wrong", "DLR host != published"), ("DLR recovery", "DLR = published"), (" truth", " published hosts"),
         ("Host recovery", "Agreement with the published host"), ("(true)", "(published)"), ("DP2", "DES-SN5YR"),
         ("real SNe", "DES SNe"), ("(no truth;", "(no reference;")]


def relabel_figures():
    """Rewrite the truth wording of every figure when it is saved."""
    import matplotlib.figure
    import matplotlib.text
    import matplotlib.ticker
    save = matplotlib.figure.Figure.savefig

    def sub(s):
        for old, new in WORDS:
            s = s.replace(old, new)
        return s

    def savefig(fig, *a, **k):
        for t in fig.findobj(matplotlib.text.Text):
            t.set_text(sub(t.get_text()))
        for ax in fig.axes:   # fixed tick labels are redrawn from their formatter
            for axis in (ax.xaxis, ax.yaxis):
                fmt = axis.get_major_formatter()
                if isinstance(fmt, (matplotlib.ticker.FixedFormatter, matplotlib.ticker.FuncFormatter)):
                    axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos, f=fmt: sub(f(x, pos))))
        return save(fig, *a, **k)
    matplotlib.figure.Figure.savefig = savefig


def frac(mask, n=None):
    k = int(np.sum(mask)); n = len(mask) if n is None else n
    p, lo, hi = wilson(k, n)
    return {"k": k, "n": n, "frac": float(p), "ci68": [float(lo), float(hi)]}


def recovery(a):
    """des_recovery: dp02_fig.recovery (agreement with the published host in the truth framing) plus HOSTGAL2,
    the published host's DR2 class, per-band / per-profile agreement from the sweep and the published d_DLR."""
    import dp02_fig
    df, rnd = dp02_fig.load(a.associations), dp02_fig.load(a.random)
    with tempfile.TemporaryDirectory() as d:
        dp02_fig.recovery(df, f"{d}/base.json", a.hostless_cut, rnd, a.p_threshold)
        with open(f"{d}/base.json") as f:
            res = json.load(f)
    cut, t = HOSTLESS_CUT[a.hostless_cut], P_TRUE[a.p_threshold]
    s = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    df = df.merge(s.drop(columns=[c for c in s.columns if c in df.columns and c != "sn_id"]), on="sn_id", how="left")
    full = df
    st_full = full.true_host_status.to_numpy()
    df = df[L.rung_mask(df, 4)]   # every agreement number below on the headline's denominator (rung 4), 2 Nov 2026
    n = len(df)
    sb, dl = df.sbw_objectId.to_numpy(), df.dlr_objectId.to_numpy()
    h1, h2 = df.host_objectId.to_numpy(), df.host2_objectId.to_numpy()
    hl = sbw_hostless(df, cut)
    ok = sbw_success(df, t)
    rows, cols = prob_rows(df, cut, t), dlr_cols(df)
    res["reference"] = ("published DES-SN5YR host (HOSTGAL), the DES-SN DLR assignment on deep SN-field coadds, "
                        "matched to DES DR2 under reference_match: agreement with a DLR-derived reference, not recovery of truth")
    res["agreement_with_published"] = {
        "sbw": frac(ok), "sbw_hosted": frac(ok & ~hl),
        "sbw_note": f"P_SBW(published host) >= {t:g}: the published host is in SBW's credible set",
        "dlr": frac(df.dlr_correct), "mean_p_sbw_published": res["mean_p_sbw_true"],
        "published_matched_in_dr2": frac(h1 >= 0), "published_kept_as_candidate": frac(df.true_in_cand)}
    res["p_band_matrix_fractions"] = {f"{r}/{c}": float((rm & cm).mean()) for r, rm in rows.items() for c, cm in cols.items()}
    res["p_band_matrix_ci68"] = {f"{r}/{c}": frac(rm & cm)["ci68"] for r, rm in rows.items() for c, cm in cols.items()}
    res["p_bands_ci68"] = {r: frac(rm)["ci68"] for r, rm in rows.items()}
    has2 = h2 >= 0
    res["hostgal2"] = {
        "n_present": int(df.host2_ra.notna().sum()), "n_matched": int(has2.sum()),
        "sbw_argmax_is_host2": frac((sb == h2) & has2), "dlr_is_host2": frac((dl == h2) & has2),
        "sbw_argmax_is_host_or_host2": frac(((sb == h1) & (h1 >= 0)) | ((sb == h2) & has2)),
        "dlr_is_host_or_host2": frac(((dl == h1) & (h1 >= 0)) | ((dl == h2) & has2)),
        "of_sbw_ne_published_hosted_is_host2": frac(((sb == h2) & has2)[(sb >= 0) & (sb != h1)]),
        "of_dlr_ne_published_hosted_is_host2": frac(((dl == h2) & has2)[(dl >= 0) & (dl != h1)])}
    st = df.true_host_status.to_numpy()
    res["published_host_status"] = {k: int((st_full == k).sum()) for k in L.STATUS}   # the ladder, over the whole sample
    res["published_host_status_note"] = "over all sample SNe; every agreement number is on rung 4 (published host usable)"
    ecc = full.host_extended_class_coadd.to_numpy()
    res["published_host_extended_class_coadd"] = {str(k): int((ecc == k).sum()) for k in sorted(set(ecc.tolist()))}
    dd = df.host_ddlr_published.to_numpy(float)
    res["published_ddlr"] = {"median": float(np.nanmedian(dd)), "frac_gt_3": float(np.mean(dd > 3)),
                             "median_when_dlr_agrees": float(np.nanmedian(dd[df.dlr_correct])),
                             "median_when_dlr_disagrees": float(np.nanmedian(dd[~df.dlr_correct])) if (~df.dlr_correct).any() else None,
                             "our_d_dlr_of_published_median": float(np.nanmedian(df.d_dlr_true))}
    res["sbw_ne_published_split"] = {
        "note": "the first three over the whole sample, the kept_* over rung 4",
        "published_absent_from_dr2_or_search": int((st_full == "absent").sum()),
        "published_point_source": int((st_full == "point_source").sum()),
        "published_unusable": int((st_full == "unusable").sum()),
        "kept_sbw_hostless": int(((st == "kept") & hl).sum()),
        "kept_sbw_excludes_published": int(((st == "kept") & ~hl & ~ok).sum()),
        "kept_sbw_excludes_published_and_dlr_elsewhere": int(((st == "kept") & ~hl & ~ok & (dl >= 0) & (dl != h1)).sum())}
    sw = pq.read_table(a.sweep).to_pandas(ignore_metadata=True).merge(s[["sn_id", "host_objectId"]], on="sn_id")
    sw = sw[sw.sn_id.isin(df.sn_id[L.rung_mask(df, 4)])]   # the headline's denominator (rung 4), 2 Nov 2026
    res["per_band_profile"] = {"note": "on the headline's denominator: SNe whose published host is usable by both methods (rung 4)"}
    for (b, p), g in sw.groupby(["band", "profile"]):
        hb = sbw_hostless(g, cut)
        ob = sbw_success(g, t)
        res["per_band_profile"][f"{b}/{p}"] = {
            "sbw_agrees": float(ob.mean()), "sbw_hosted_agrees": float((ob & ~hb).mean()),
            "mean_p_sbw_published": float(L.p_true(g).mean()),
            "dlr_agrees": float((g.dlr_objectId == g.host_objectId).mean()),
            "sbw_hostless": float(hb.mean()), "dlr_hostless": float((g.dlr_objectId < 0).mean())}
    with open(a.output, "w") as f:
        json.dump(res, f, indent=2, default=float)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    kind, argv = sys.argv[1], sys.argv[2:]
    if kind == "associate":
        import dp02_associate
        dp02_associate.main(argv)
    elif kind == "sweep":
        import dp02_sweep
        dp02_sweep.BANDS = SWEEP_BANDS
        dp02_sweep.main(argv)
    elif kind == "fig":
        import dp02_fig
        dp02_fig.TIER["des"] = "DES-SN5YR"
        relabel_figures()
        dp02_fig.main(argv)
    elif kind == "dp2fig":
        import dp2_fig
        dp2_fig.BANDS = SWEEP_BANDS
        dp2_fig.PROFILE_SHORT.update({"exponential": "exponential (FLUX_RADIUS)", "sersic": "n = 4 (FLUX_RADIUS)",
                                      "moments": "moment Gaussian (A, B)"})
        relabel_figures()
        dp2_fig.main(argv)
    elif kind == "recovery":
        ap = argparse.ArgumentParser(description=recovery.__doc__)
        for k in ("associations", "random", "sample", "sweep", "output"):
            ap.add_argument(f"--{k}", required=True)
        ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
        ap.add_argument("--p-threshold", required=True, choices=sorted(P_TRUE))
        recovery(ap.parse_args(argv))
    else:
        raise SystemExit(f"unknown kind {kind!r}\n{__doc__}")


if __name__ == "__main__":
    main()
