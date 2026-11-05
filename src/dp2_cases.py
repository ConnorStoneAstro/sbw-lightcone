"""DP2 disagreement cases for the cutout gallery: a seeded draw of transients per class, with the positions
and half-light ellipses of the SBW most probable host and the DLR host, and in sbw_set the other members of SBW's
credible set (P_SBW >= t) with their P and half-light ellipses.

Classes (headline band and profile; SBW hostless = flb_best < cut, DLR hostless = no d_DLR <= dlr_max):
sbw_ne_dlr (both name a host, different objects), sbw_hostless_only, dlr_hostless_only, and band_flip (SBW
names a host in the headline band and in another band, and they differ; the other band's host is drawn too,
from the dp2_association_sweep of the same profile). Ellipses are in arcsec with the angle counter-clockwise
from +x = West toward North (the sky_to_plane frame): the SBW ellipse is the profile's half-light ellipse, the
DLR ellipse is the one d_DLR is measured on (d_DLR = 1 on it).
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import dp02_lib as L
import dp2_lib
from dp02_lib import BANDS, HOSTLESS_CUT, P_TRUE, sbw_credible_sets, sbw_hostless

N_PER_CLASS = 12
SEED = 20261028


def ellipses(cand, band, profile, dlr_radius):
    """Per candidate row: SBW half-light ellipse (a, b, phi) and DLR ellipse (a, b, phi)."""
    a, _, _, q, phi, _ = L.light_model(cand, band, profile)
    g = lambda k: cand[f"{band}_{k}"].to_numpy(dtype=float)
    qm, _, l1, l2, _ = L.moment_shape(g("ixx"), g("iyy"), g("ixy"), g("ixxPSF"), g("iyyPSF"), g("ixyPSF"))
    da, db, _, dphi = L.dlr_ellipse(cand, band, dlr_radius, qm, l1, l2)
    return pd.DataFrame({"sbw_a": a, "sbw_b": a * q, "sbw_phi": phi, "dlr_a": da, "dlr_b": db, "dlr_phi": dphi},
                        index=cand.index)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--exclude-suspect", choices=["keep", "drop"], default="keep")
    ap.add_argument("--p-threshold", required=True, choices=sorted(P_TRUE), help="credible set: P_SBW >= t")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    dp2_lib.install(a.exclude_suspect)
    t = pq.read_table(a.associations)
    df = t.to_pandas(ignore_metadata=True)
    attrs = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if not k.startswith(b"pandas")}
    band, profile, dlr_radius = attrs["band"], attrs["profile"], attrs["dlr_radius"]
    cut = HOSTLESS_CUT[a.hostless_cut]
    sh, dh = sbw_hostless(df, cut), df.dlr_objectId.to_numpy() < 0
    s, d = df.sbw_objectId.to_numpy(), df.dlr_objectId.to_numpy()
    # band flips: SBW hosted in the headline band and in band b with another host; the first such band in BANDS order
    sw = pq.read_table(a.sweep).to_pandas(ignore_metadata=True)
    sw = sw[sw.profile == profile].set_index(["band", "sn_id"])
    flip_band = pd.Series("", index=df.sn_id)
    flip_id = pd.Series(np.int64(-1), index=df.sn_id)
    for b in BANDS:
        if b == band:
            continue
        x = sw.loc[b].reindex(df.sn_id)
        f = ~sh & ~sbw_hostless(x, cut) & (x.sbw_objectId.to_numpy() != s) & (flip_band.to_numpy() == "")
        flip_band[f] = b
        flip_id[f] = x.sbw_objectId.to_numpy()[f]
    df["flip_band"], df["flip_objectId"] = flip_band.to_numpy(), flip_id.to_numpy()
    classes = {"sbw_ne_dlr": ~sh & ~dh & (s != d), "sbw_hostless_only": sh & ~dh, "dlr_hostless_only": ~sh & dh,
               "band_flip": df.flip_band.to_numpy() != ""}
    rng = np.random.default_rng(SEED)
    picks = []
    for name, sel in classes.items():
        idx = np.flatnonzero(sel)
        idx = rng.choice(idx, size=min(N_PER_CLASS, len(idx)), replace=False)
        picks.append(df.iloc[np.sort(idx)].assign(**{"class": name, "class_n": int(sel.sum())}))
    c = pd.concat(picks, ignore_index=True)

    cand = L.load_candidates(a.candidates, bands=band)
    cand = cand[cand.sn_id.isin(c.sn_id)].drop_duplicates(["sn_id", "objectId"]).reset_index(drop=True)
    thr = P_TRUE[a.p_threshold]
    c["sbw_set"] = sbw_credible_sets(c.assign(host_objectId=c.get("host_objectId", -1)), cand, attrs, thr,
                                     with_ellipse=True)
    c["sbw_set_t"] = thr
    cand = pd.concat([cand[["sn_id", "objectId", "coord_ra", "coord_dec"]], ellipses(cand, band, profile, dlr_radius)],
                     axis=1).set_index(["sn_id", "objectId"])
    for name, col in [("sbw", "sbw_objectId"), ("dlr", "dlr_objectId"), ("flip", "flip_objectId")]:
        pos = cand.reindex(pd.MultiIndex.from_arrays([c.sn_id.to_numpy(), c[col].to_numpy()]))
        c[f"{name}_ra"], c[f"{name}_dec"] = pos.coord_ra.to_numpy(), pos.coord_dec.to_numpy()
        kind = "dlr" if name == "dlr" else "sbw"
        for k in ("a", "b", "phi"):
            c[f"{name}_{k}"] = pos[f"{kind}_{k}"].to_numpy()
    c["case_id"] = [f"{k}_{i:02d}" for k, i in zip(c["class"], c.groupby("class").cumcount())]
    c["label"] = [f"{k} | {sid} | P(SBW) {p:.2f}, P(DLR host) {pd_:.2f} | flb {f:.2g} | d_DLR {dd:.1f}"
                  + (f" | {fb} host differs" if fb else "")
                  for k, sid, p, pd_, f, dd, fb in zip(c["class"], c.sn_id, c.p_sbw_max, c.p_sbw_dlr, c.flb_best,
                                                       c.d_dlr_min, c.flip_band)]
    cols = (["case_id", "class", "class_n", "sn_id", "sn_ra", "sn_dec", "label"]
            + [f"{n}_{k}" for n in ("sbw", "dlr", "flip") for k in ("ra", "dec", "a", "b", "phi")]
            + ["sbw_set", "sbw_set_t", "p_sbw_max", "p_sbw_dlr", "flb_best", "d_dlr_min", "d_dlr_sbw", "sep_sbw", "sep_dlr", "n_cand_used",
               "flip_band", "sbw_objectId", "dlr_objectId", "flip_objectId", "z"])
    c[[k for k in cols if k in c.columns]].to_csv(a.output, index=False, float_format="%.10g")


if __name__ == "__main__":
    main()
