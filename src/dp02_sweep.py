"""Sweeps that recompute the association over bands or profiles:
calibration-by-band, calibration-by-profile, band dependence, profile table."""
from __future__ import annotations

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from dp02_fig import plot_calibration, plt as _plt  # noqa: F401  (rcParams)
from dp02_lib import BANDS, COLORS, DLR_MAX, PROFILES, PROFILE_LABEL, QUALITY, RADIUS, associate_all, p_true, calibration_curve, load_candidates, wilson


def run(sample, candpath, bands, profiles, cfg, workers):
    res = {}
    for b in bands:
        cand = load_candidates(candpath, bands=b)
        for p in profiles:
            res[(b, p)] = associate_all(sample, cand, band=b, profile=p, workers=workers, **cfg)
        del cand
    return res


def panel_grid(res, keys, titles, out, suptitle):
    n = len(keys)
    ncol = 3 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 3.0 * nrow), squeeze=False)
    for ax, k, t in zip(axes.flat, keys, titles):
        df = res[k]
        _, ece = plot_calibration(ax, df)
        ax.set_title(f"{t}\nmean P(true) {100 * p_true(df).mean():.2f}%, most probable = true "
                     f"{100 * df.sbw_correct.mean():.2f}%  N={len(df):,}", fontsize=7)
        ax.set_xlabel(""); ax.set_ylabel("")
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.supxlabel("$P_{\\mathrm{SBW}}$ of the most probable host", fontsize=9)
    fig.supylabel("fraction where the most probable host is the true host", fontsize=9)
    fig.suptitle(suptitle, fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def band_dependence(res, profile, out, cfgtxt):
    bands = list(BANDS)
    rec = [res[(b, profile)].sbw_correct.mean() for b in bands]
    mp = [p_true(res[(b, profile)]).mean() for b in bands]
    drec = [res[(b, profile)].dlr_correct.mean() for b in bands]
    agr = [res[(b, profile)].agree.mean() for b in bands]
    n = [len(res[(b, profile)]) for b in bands]
    ref = res[("r", profile)].set_index("sn_id")["sbw_objectId"]
    changed = []
    for b in bands:
        s = res[(b, profile)].set_index("sn_id")["sbw_objectId"]
        both = (ref.reindex(s.index) >= 0) & (s >= 0)
        changed.append(float((ref.reindex(s.index)[both] != s[both]).mean()))
    x = np.arange(6)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for y, name, color, mk, ls in [(mp, "mean $P_{\\mathrm{SBW}}$(true host)", COLORS["sbw"], "o", "-"),
                                   (rec, "SBW most probable host = true host", COLORS["sbw"], "o", "--"),
                                   (drec, "DLR host = true host", COLORS["dlr"], "s", "-"),
                                   (agr, "SBW most probable host = DLR host", COLORS["agree"], "^", "-"),
                                   (changed, "SBW most probable host differs from r band", "#7b4ea3", "d", "-")]:
        ax.plot(x, y, marker=mk, ms=5, color=color, ls=ls, label=name, mfc="none" if ls == "--" else color)
    ax.set_xticks(x, bands); ax.set_xlabel("band"); ax.set_ylabel("fraction"); ax.set_ylim(0, 1.02)
    ax.legend(fontsize=7, loc="center right")
    ax.set_title(f"Band dependence ({PROFILE_LABEL[profile]})\n{cfgtxt}", fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def profile_table(res, band, out):
    rows = []
    for p in PROFILES:
        df = res[(band, p)]
        n = len(df)
        _, ece = calibration_curve(df)
        ps, lo, hi = wilson(df.sbw_correct.sum(), n)
        rows.append({"profile": p, "label": PROFILE_LABEL[p], "band": band, "n": n,
                     "mean_p_sbw_true": p_true(df).mean(),
                     "sbw_recovery": ps, "sbw_recovery_lo": lo, "sbw_recovery_hi": hi,
                     "dlr_recovery": df.dlr_correct.mean(), "agreement": df.agree.mean(),
                     "sbw_no_host": (df.sbw_objectId < 0).mean(), "median_p_true": np.nanmedian(df.p_sbw_true),
                     "expected_calibration_error": ece,
                     "disagree_sbw_right": (~df.agree & df.sbw_correct).mean(),
                     "disagree_dlr_right": (~df.agree & df.dlr_correct).mean()})
    pd.DataFrame(rows).to_csv(out, index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["band", "profile", "band_dependence", "profile_table"])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--profile", choices=PROFILES + ("stamps",))
    ap.add_argument("--band", choices=list(BANDS))
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--dlr-max", required=True, choices=sorted(DLR_MAX))
    ap.add_argument("--dlr-radius", required=True, choices=["half_light", "moments"])
    ap.add_argument("--extendedness", required=True, choices=["ext05", "none"])
    ap.add_argument("--candidate-quality", required=True, choices=QUALITY)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--max-sn", type=int, default=None, help="probe: subsample the SNe")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    cfg = dict(search_radius=a.search_radius, dlr_max=a.dlr_max, dlr_radius=a.dlr_radius, extendedness=a.extendedness,
               quality=a.candidate_quality)
    cfgtxt = f"{a.search_radius} | d_DLR<={a.dlr_max[1:]} ({a.dlr_radius}) | {a.extendedness} | quality {a.candidate_quality}"
    sample = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    if a.max_sn:
        sample = sample.sample(a.max_sn, random_state=0)
    if a.mode in ("band", "band_dependence"):
        if a.profile == "stamps":
            raise SystemExit("stamps not available on DP0.2")
        res = run(sample, a.candidates, BANDS, [a.profile], cfg, a.workers)
        if a.mode == "band":
            panel_grid(res, [(b, a.profile) for b in BANDS], [f"{b} band" for b in BANDS], a.output,
                       f"SBW calibration per band ({PROFILE_LABEL[a.profile]}) | {cfgtxt}")
        else:
            band_dependence(res, a.profile, a.output, cfgtxt)
    else:
        res = run(sample, a.candidates, [a.band], PROFILES, cfg, a.workers)
        if a.mode == "profile":
            panel_grid(res, [(a.band, p) for p in PROFILES], [PROFILE_LABEL[p] for p in PROFILES], a.output,
                       f"SBW calibration per light model ({a.band} band) | {cfgtxt}")
        else:
            profile_table(res, a.band, a.output)


if __name__ == "__main__":
    main()
