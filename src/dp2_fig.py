"""DP2 (no truth) summaries: SBW-DLR disagreement in the P_SBW-band framing, agreement versus d_DLR, band and
profile consistency, and the top-P_SBW / credible-set distributions against DP0.2.

Without a true host, SBW is scored against DLR's pick: P_SBW of the DLR host falls in a band (> 0.5,
t..0.5, < t), or SBW calls the transient hostless (flb_best < cut); DLR either names a host (the same object
as SBW's most probable host or another) or none (no candidate with d_DLR <= dlr_max).
"""
from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from dp02_fig import DP02_OPTIMA, binned_fraction, binned_mean, load, set_optima, subtitle
from dp02_lib import BANDS, COLORS, HOSTLESS_CUT, P_TRUE, PROFILES, sbw_hostless

PROFILE_SHORT = {"exponential": "exponential fit", "sersic": "free-n Sersic fit", "moments": "moment Gaussian"}
CATS = ("same_host", "different_host", "sbw_hostless_only", "dlr_hostless_only", "both_hostless")
CAT_LABEL = {"same_host": "same host", "different_host": "different hosts", "sbw_hostless_only": "SBW: no host only",
             "dlr_hostless_only": "DLR: no host only", "both_hostless": "both: no host"}
CAT_COLOR = {"same_host": COLORS["agree"], "different_host": "#7b4ea3", "sbw_hostless_only": COLORS["sbw"],
             "dlr_hostless_only": COLORS["dlr"], "both_hostless": "#bbbbbb"}


def categories(df, cut):
    """SBW (hostless at flb_best < cut) versus DLR (hostless beyond dlr_max) outcome per transient."""
    sh = sbw_hostless(df, cut)
    dh = df["dlr_objectId"].to_numpy() < 0
    same = df["sbw_objectId"].to_numpy() == df["dlr_objectId"].to_numpy()
    return {"same_host": ~sh & ~dh & same, "different_host": ~sh & ~dh & ~same,
            "sbw_hostless_only": sh & ~dh, "dlr_hostless_only": ~sh & dh, "both_hostless": sh & dh}


def p_band_cells(df, cut, t):
    """Rows: P_SBW of the DLR host (> 0.5 / t..0.5 / < t) or SBW hostless; columns: DLR host is / is not SBW's
    most probable host, or DLR has none (rows then just SBW host / SBW hostless)."""
    sh = sbw_hostless(df, cut)
    dh = df["dlr_objectId"].to_numpy() < 0
    same = df["sbw_objectId"].to_numpy() == df["dlr_objectId"].to_numpy()
    p = np.nan_to_num(df["p_sbw_dlr"].to_numpy(dtype=float), nan=0.0)
    rows = {"p_dlr_gt_0p5": ~sh & (p > 0.5), "p_dlr_t_to_0p5": ~sh & (p >= t) & (p <= 0.5),
            "p_dlr_lt_t": ~sh & (p < t), "sbw_hostless": sh}
    cols = {"dlr_eq_sbw": ~dh & same, "dlr_ne_sbw": ~dh & ~same}
    out = {f"{r}/{c}": int((rm & cm).sum()) for r, rm in rows.items() for c, cm in cols.items()}
    out["sbw_host/dlr_hostless"] = int((~sh & dh).sum())
    out["sbw_hostless/dlr_hostless"] = int((sh & dh).sum())
    out.update(n=len(df), p_threshold=t, sbw_hostless_cut=cut)
    return out


def summary(df, cut, t):
    n = len(df)
    c = {k: float(v.mean()) for k, v in categories(df, cut).items()}
    p = df["p_sbw_dlr"].to_numpy(dtype=float)
    hosted = df["dlr_objectId"].to_numpy() >= 0
    return {"n": n, **c,
            "disagree_any": c["different_host"] + c["sbw_hostless_only"] + c["dlr_hostless_only"],
            "argmax_disagree": float((df.sbw_objectId.to_numpy() != df.dlr_objectId.to_numpy()).mean()),
            "sbw_hostless": float(sbw_hostless(df, cut).mean()), "dlr_hostless": float((~hosted).mean()),
            "mean_p_sbw_of_dlr_host": float(np.nanmean(p[hosted])) if hosted.any() else np.nan,
            "p_sbw_max_median_hosted": float(np.nanmedian(df["p_sbw_max"].to_numpy(dtype=float)[~sbw_hostless(df, cut)])) if (~sbw_hostless(df, cut)).any() else np.nan,
            "p_sbw_max_gt_0p9_hosted": float((df["p_sbw_max"].to_numpy(dtype=float)[~sbw_hostless(df, cut)] > 0.9).mean()) if (~sbw_hostless(df, cut)).any() else np.nan,
            "frac_p_sbw_of_dlr_host_lt_t": float((np.nan_to_num(p[hosted]) < t).mean()) if hosted.any() else np.nan}


def read_sweep(path):
    t = pq.read_table(path)
    return t.to_pandas(ignore_metadata=True)


def disagreement_rate(df, sweep, rnd, out, hostless_cut, p_threshold):
    cut, t = HOSTLESS_CUT[hostless_cut], P_TRUE[p_threshold]
    res = {"definition": __doc__.split("\n\n")[1].replace("\n", " "),
           "headline": summary(df, cut, t), "p_band_matrix": p_band_cells(df, cut, t),
           "random_positions": {**summary(rnd, cut, t), "p_band_matrix": p_band_cells(rnd, cut, t)},
           "per_band": {}, "per_profile": {}, "band_flip": {}}
    prof, band = df.attrs.get("profile"), df.attrs.get("band")
    for b in BANDS:
        res["per_band"][b] = summary(sweep[(sweep.band == b) & (sweep.profile == prof)], cut, t)
    for p in PROFILES:
        res["per_profile"][p] = summary(sweep[(sweep.band == band) & (sweep.profile == p)], cut, t)
    for p in PROFILES:
        res["band_flip"][p] = band_changes(sweep, p, band, cut)
    res["profile_agreement"] = profile_agreement(sweep, band, cut)
    res["settings"] = df.attrs
    with open(out, "w") as f:
        json.dump(res, f, indent=2, default=float)


def band_changes(sweep, profile, ref, cut):
    """Per band: fraction of transients hosted by SBW in both bands whose most probable host differs from the
    reference band's, the same for DLR, and the fraction that changes hostless status."""
    s = sweep[sweep.profile == profile].set_index(["band", "sn_id"])
    r = s.loc[ref]
    out = {}
    for b in BANDS:
        x = s.loc[b].reindex(r.index)
        hs_r, hs_b = sbw_hostless(r, cut), sbw_hostless(x, cut)
        both = ~hs_r & ~hs_b
        dboth = (r.dlr_objectId.to_numpy() >= 0) & (x.dlr_objectId.to_numpy() >= 0)
        out[b] = {"sbw_host_changes": float((r.sbw_objectId.to_numpy() != x.sbw_objectId.to_numpy())[both].mean()),
                  "dlr_host_changes": float((r.dlr_objectId.to_numpy() != x.dlr_objectId.to_numpy())[dboth].mean()),
                  "sbw_hostless_status_changes": float((hs_r != hs_b).mean()),
                  "dlr_hostless_status_changes": float(((r.dlr_objectId.to_numpy() < 0) != (x.dlr_objectId.to_numpy() < 0)).mean()),
                  "n_both_hosted": int(both.sum())}
    return out


def profile_agreement(sweep, band, cut):
    s = sweep[sweep.band == band].set_index(["profile", "sn_id"])
    out = {}
    for i, p in enumerate(PROFILES):
        for q in PROFILES[i + 1:]:
            a, b = s.loc[p], s.loc[q].reindex(s.loc[p].index)
            both = ~sbw_hostless(a, cut) & ~sbw_hostless(b, cut)
            out[f"{p}_vs_{q}"] = {"same_sbw_host": float((a.sbw_objectId.to_numpy() == b.sbw_objectId.to_numpy())[both].mean()),
                                  "same_hostless_status": float((sbw_hostless(a, cut) == sbw_hostless(b, cut)).mean()),
                                  "n_both_hosted": int(both.sum())}
    return out


def agreement_vs_ddlr(df, out, hostless_cut):
    cut = HOSTLESS_CUT[hostless_cut]
    agree = df.agree.to_numpy()
    fig, axes = plt.subplots(2, 3, figsize=(11, 4.8), height_ratios=[3, 1], sharex="col")
    for j, (col, lab, edges, logx) in enumerate([
            ("d_dlr_min", "$d_{DLR}$ of the DLR host (smallest $d_{DLR}$)", np.linspace(0, 8, 17), False),
            ("d_dlr_sbw", "$d_{DLR}$ of the SBW most probable host", np.linspace(0, 8, 17), False),
            ("p_sbw_max", "$P_{\\mathrm{SBW}}$ of the most probable host", np.array([0, .3, .5, .6, .7, .8, .9, .95, .99, 1.0001]), False)]):
        x = df[col].to_numpy(dtype=float)
        ax = axes[0, j]
        c, p, lo, hi, n = binned_fraction(x, agree, edges)
        ax.fill_between(c, lo, hi, color=COLORS["agree"], alpha=0.2)
        ax.plot(c, p, "o-", ms=3, color=COLORS["agree"], label="SBW most probable host = DLR host")
        c2, m2, _ = binned_mean(x, np.nan_to_num(df.p_sbw_dlr.to_numpy(dtype=float), nan=0.0), edges)
        ax.plot(c2, m2, "s-", ms=3, color=COLORS["dlr"], label="mean $P_{\\mathrm{SBW}}$ of the DLR host (0 if none)")
        c3, p3, _, _, _ = binned_fraction(x, sbw_hostless(df, cut), edges)
        ax.plot(c3, p3, "^-", ms=3, color=COLORS["sbw"], label=f"SBW: no host (flb < {cut:g})")
        ax.set_ylim(0, 1.02); ax.set_ylabel("fraction")
        if j < 2:
            ax.axvline(float(df.attrs.get("dlr_max", "d4")[1:]), color=COLORS["ref"], ls=":", lw=1)
        if j == 0:
            ax.legend(fontsize=7, loc="lower left")
        axes[1, j].bar(c, n, width=np.diff(edges), color=COLORS["ref"], alpha=0.5)
        axes[1, j].set_xlabel(lab); axes[1, j].set_ylabel("transients")
    fig.suptitle(f"DP2: SBW-DLR agreement versus $d_{{DLR}}$ and $P_{{\\mathrm{{SBW}}}}$ (no truth; overall agreement "
                 f"{100 * agree.mean():.1f}%)\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def band_consistency(df, sweep, out, hostless_cut):
    cut = HOSTLESS_CUT[hostless_cut]
    prof, ref = df.attrs.get("profile"), df.attrs.get("band")
    x = np.arange(len(BANDS))
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.0))
    (aa, ab), (ac, ad) = axes
    for p, ls in zip(PROFILES, ("-", "--", ":")):
        bc = band_changes(sweep, p, ref, cut)
        aa.plot(x, [bc[b]["sbw_host_changes"] for b in BANDS], "o", ls=ls, color=COLORS["sbw"], label=f"SBW, {PROFILE_SHORT[p]}")
    bc = band_changes(sweep, prof, ref, cut)
    aa.plot(x, [bc[b]["dlr_host_changes"] for b in BANDS], "s-", color=COLORS["dlr"], label="DLR")
    aa.set_xticks(x, BANDS); aa.set_ylabel(f"host differs from the {ref} band"); aa.set_ylim(bottom=0)
    aa.legend(fontsize=7); aa.set_title(f"(a) host change versus {ref} band (both hosted)", fontsize=9, loc="left")
    # (b) categories per band, headline profile
    bottom = np.zeros(len(BANDS))
    for k in CATS:
        v = np.array([categories(sweep[(sweep.band == b) & (sweep.profile == prof)], cut)[k].mean() for b in BANDS])
        ab.bar(x, v, bottom=bottom, color=CAT_COLOR[k], label=CAT_LABEL[k]); bottom += v
    ab.set_xticks(x, BANDS); ab.set_ylim(0, 1); ab.set_ylabel("fraction of transients")
    ab.legend(fontsize=7, loc="lower right"); ab.set_title(f"(b) SBW vs DLR per band ({PROFILE_SHORT[prof]})", fontsize=9, loc="left")
    # (c) hostless fractions per band and profile
    for p, ls in zip(PROFILES, ("-", "--", ":")):
        aa_ = [sbw_hostless(sweep[(sweep.band == b) & (sweep.profile == p)], cut).mean() for b in BANDS]
        ac.plot(x, aa_, "o", ls=ls, color=COLORS["sbw"], label=f"SBW no host, {PROFILE_SHORT[p]}")
    ac.plot(x, [(sweep[(sweep.band == b) & (sweep.profile == prof)].dlr_objectId < 0).mean() for b in BANDS], "s-",
            color=COLORS["dlr"], label="DLR no host")
    ac.set_xticks(x, BANDS); ac.set_ylim(bottom=0); ac.set_ylabel("fraction hostless")
    ac.legend(fontsize=7); ac.set_title(f"(c) hostless calls (flb < {cut:g}; $d_{{DLR}}$ > {df.attrs.get('dlr_max', 'd4')[1:]})", fontsize=9, loc="left")
    # (d) profile agreement per band
    for (p, q), ls in zip([("exponential", "sersic"), ("exponential", "moments"), ("sersic", "moments")], ("-", "--", ":")):
        v = [profile_agreement(sweep, b, cut)[f"{p}_vs_{q}"]["same_sbw_host"] for b in BANDS]
        ad.plot(x, v, "o", ls=ls, color=COLORS["agree"], label=f"{PROFILE_SHORT[p]} vs {PROFILE_SHORT[q]}")
    ad.set_xticks(x, BANDS); ad.set_ylabel("same SBW most probable host"); ad.set_ylim(top=1.005)
    ad.legend(fontsize=7, loc="lower center"); ad.set_title("(d) SBW host agreement between light models (both hosted)", fontsize=9, loc="left")
    fig.suptitle(f"DP2 band and light-model consistency, {sweep.sn_id.nunique():,} transients\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def p_top_matrix(df):
    tops = df["p_top"].to_numpy()
    L = max(1, max(len(x) for x in tops))
    P = np.zeros((len(tops), L), np.float32)
    for i, x in enumerate(tops):
        P[i, :len(x)] = x
    return P


T_GRID = np.geomspace(1e-3, 0.5, 61)


def p_top_host(df, rnd, ref, out, hostless_cut):
    cut = HOSTLESS_CUT[hostless_cut]
    sets = [("DP2 transients", df, COLORS["sbw"], "-"), ("DP2 random positions", rnd, COLORS["agree"], "-"),
            ("DP0.2 SNe (reference)", ref, "black", "--")]
    fig, (aa, ab) = plt.subplots(1, 2, figsize=(9.6, 3.8))
    edges = np.linspace(0, 1, 41)
    for lab, d, col, ls in sets:
        keep = ~sbw_hostless(d, cut)
        p = d.p_sbw_max.to_numpy(dtype=float)[keep]
        aa.hist(p, bins=edges, density=True, histtype="step", color=col, ls=ls, lw=1.5,
                label=f"{lab}: N = {keep.sum():,}, median {np.median(p):.3f}, P > 0.9 {100 * (p > 0.9).mean():.1f}%")
        P = p_top_matrix(d[keep])
        ab.plot(T_GRID, [(P > t).sum(1).mean() for t in T_GRID], color=col, ls=ls, label=lab)
    aa.set_yscale("log"); aa.set_xlabel("$P_{\\mathrm{SBW}}$ of the most probable host"); aa.set_ylabel("density")
    aa.legend(fontsize=6.5, loc="upper left"); aa.set_title(f"(a) top $P_{{\\mathrm{{SBW}}}}$, SBW-hosted (flb $\\geq$ {cut:g})", fontsize=9, loc="left")
    ab.set_xscale("log"); ab.axvline(DP02_OPTIMA["credible_set_t"], color=COLORS["ref"], ls=":", lw=1)
    ab.text(DP02_OPTIMA["credible_set_t"] * 1.08, ab.get_ylim()[1] * 0.97, f"DP0.2 knee t = {DP02_OPTIMA['credible_set_t']:g}",
            fontsize=7, va="top")
    ab.set_xlabel("threshold t (set = candidates with $P_{\\mathrm{SBW}}$ > t)"); ab.set_ylabel("mean credible-set size")
    ab.legend(fontsize=7); ab.set_title("(b) credible-set size versus threshold", fontsize=9, loc="left")
    fig.suptitle("DP2 SBW probabilities without truth, DP0.2 dashed for reference\n" + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def agreement_matrix(df, rnd, out, hostless_cut, p_threshold):
    """Figure form of the P-band matrix (dp2_disagreement_rate) in the layout of dp02_agreement_matrix.
    (a) rows: P_SBW of the DLR host > 0.5 / t..0.5 / < t, or SBW hostless (flb_best < cut); columns: the DLR host is
    SBW's most probable host / another object / DLR hostless (where the P-band rows do not apply, the DLR-hostless
    column holds SBW hosted and SBW hostless). (b) SBW x DLR hosted / hostless for the transients where SBW's most
    probable host has P_SBW > t. (c) the same 2x2 on random positions, where hostless is the correct answer."""
    cut, t = HOSTLESS_CUT[hostless_cut], P_TRUE[p_threshold]
    n = len(df)
    sh = sbw_hostless(df, cut)
    d = df["dlr_objectId"].to_numpy()
    dh = d < 0
    same = df["sbw_objectId"].to_numpy() == d
    p = np.nan_to_num(df["p_sbw_dlr"].to_numpy(dtype=float), nan=0.0)
    rows = [("$P_{\\mathrm{SBW}}$(DLR host) > 0.5", ~sh & ~dh & (p > 0.5)),
            (f"$P_{{\\mathrm{{SBW}}}}$(DLR host)\nin [{t:g}, 0.5]", ~sh & ~dh & (p >= t) & (p <= 0.5)),
            (f"$P_{{\\mathrm{{SBW}}}}$(DLR host) < {t:g}", ~sh & ~dh & (p < t)),
            (f"SBW: no host\n(flb < {cut:g})", sh)]
    cols = [("DLR host = SBW most\nprobable host", ~dh & same), ("DLR host = another\nobject", ~dh & ~same), ("DLR: no host", dh)]
    M = np.array([[int((r & c).sum()) for _, c in cols] for _, r in rows], dtype=float)
    M[:3, 2] = np.nan                       # P bands undefined without a DLR host
    n_sbw_hosted_dlr_hostless = int((~sh & dh).sum())
    fig, (ax, ab, ar) = plt.subplots(1, 3, figsize=(13.2, 4.8), width_ratios=[3, 1.9, 1.9])
    ax.imshow(np.nan_to_num(M) / n, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(len(rows)):
        for j in range(len(cols)):
            if np.isnan(M[i, j]):
                continue
            ax.text(j, i, f"{int(M[i, j]):,}\n{100 * M[i, j] / n:.2f}%", ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] / n > 0.5 else "black")
    ax.add_patch(plt.Rectangle((1.5, -0.5), 1, 3, facecolor="#eeeeee", edgecolor="none", zorder=0))
    ax.text(2, 1, f"SBW assigns a host\n{n_sbw_hosted_dlr_hostless:,}\n{100 * n_sbw_hosted_dlr_hostless / n:.2f}%",
            ha="center", va="center", fontsize=8)
    ax.set_xticks(range(3), [c for c, _ in cols]); ax.set_yticks(range(len(rows)), [r for r, _ in rows])
    ax.grid(False); ax.tick_params(labelsize=7)
    ax.set_title(f"(a) DP2 transients (N = {n:,}; SBW most probable host = DLR host {100 * (same & ~sh & ~dh).mean():.2f}%)",
                 fontsize=8, loc="left")
    hosted = ~dh
    ax.set_xlabel(f"mean $P_{{\\mathrm{{SBW}}}}$(DLR host) where DLR names one {100 * np.nanmean(p[hosted]):.2f}%; "
                  f"hostless: SBW {100 * sh.mean():.2f}%, DLR {100 * dh.mean():.2f}%", fontsize=8)
    # (b) hosted / hostless 2x2 where SBW's most probable host has P > t
    conf = df["p_sbw_max"].to_numpy(dtype=float) > t
    nb = int(conf.sum())
    brows = [("SBW: no host", sh[conf]), ("SBW assigns host", ~sh[conf])]
    bcols = [("DLR: no host", dh[conf]), ("DLR assigns host", ~dh[conf])]
    B = np.array([[int((r & c).sum()) for _, c in bcols] for _, r in brows])
    ab.imshow(B / max(nb, 1), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(2):
        for j in range(2):
            ab.text(j, i, f"{B[i, j]:,}\n{100 * B[i, j] / max(nb, 1):.2f}%", ha="center", va="center", fontsize=8,
                    color="white" if B[i, j] / max(nb, 1) > 0.5 else "black")
    ab.set_xticks(range(2), [c for c, _ in bcols]); ab.set_yticks(range(2), [r for r, _ in brows]); ab.grid(False)
    ab.tick_params(labelsize=7)
    ab.set_title(f"(b) transients with max $P_{{\\mathrm{{SBW}}}}$ > {t:g} (N = {nb:,} of {n:,})", fontsize=8, loc="left")
    # (c) random positions
    hr = sbw_hostless(rnd, cut); dr = rnd["dlr_objectId"].to_numpy() < 0
    nr = len(rnd)
    R = np.array([[int((r & c).sum()) for c in (dr, ~dr)] for r in (hr, ~hr)])
    ar.imshow(R / nr, cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
    for i in range(2):
        for j in range(2):
            ar.text(j, i, f"{R[i, j]:,}\n{100 * R[i, j] / nr:.2f}%", ha="center", va="center", fontsize=8)
    ar.set_xticks(range(2), ["DLR: no host\n(correct)", "DLR assigns host\n(wrong)"])
    ar.set_yticks(range(2), ["SBW: no host\n(correct)", "SBW assigns host\n(wrong)"]); ar.grid(False)
    ar.tick_params(labelsize=7)
    ar.set_title(f"(c) random positions (N = {nr:,})", fontsize=8, loc="left")
    ar.set_xlabel(f"correctly rejected: SBW {100 * hr.mean():.2f}%, DLR {100 * dr.mean():.2f}%", fontsize=8)
    fig.suptitle("DP2 without truth: SBW probability of the DLR host, hosted / hostless agreement, and random positions\n"
                 + subtitle(df), fontsize=8)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("kind", choices=["disagreement_rate", "agreement_vs_ddlr", "band_consistency", "p_top_host", "agreement_matrix"])
    ap.add_argument("--associations", required=True, help="dp2_associations (headline band and profile)")
    ap.add_argument("--sweep", help="dp2_association_sweep (every band x profile)")
    ap.add_argument("--random", help="dp2_random_associations")
    ap.add_argument("--reference", help="dp02_associations, drawn dashed for comparison")
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--p-threshold", choices=sorted(P_TRUE))
    ap.add_argument("--dp02-credible-metrics", help="dp02_credible_set_metrics json of this universe (p_top_host: the DP0.2 knee)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    set_optima(credible_path=a.dp02_credible_metrics)
    df = load(a.associations)
    if a.kind == "disagreement_rate":
        disagreement_rate(df, read_sweep(a.sweep), load(a.random), a.output, a.hostless_cut, a.p_threshold)
    elif a.kind == "agreement_vs_ddlr":
        agreement_vs_ddlr(df, a.output, a.hostless_cut)
    elif a.kind == "band_consistency":
        band_consistency(df, read_sweep(a.sweep), a.output, a.hostless_cut)
    elif a.kind == "agreement_matrix":
        agreement_matrix(df, load(a.random), a.output, a.hostless_cut, a.p_threshold)
    elif a.kind == "p_top_host":
        ref = pq.read_table(a.reference, columns=["sbw_objectId", "p_sbw_max", "flb_best", "p_top"]).to_pandas(ignore_metadata=True)
        p_top_host(df, load(a.random), ref, a.output, a.hostless_cut)


if __name__ == "__main__":
    main()
