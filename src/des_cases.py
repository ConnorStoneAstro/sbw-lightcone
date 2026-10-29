"""DES disagreement cases for the cutout gallery: a seeded draw of SNe per class, with the published host position
and the positions and ellipses of the SBW most probable host and the DLR host.

Classes (headline band and profile; SBW hostless = flb_best < cut, DLR hostless = no d_DLR <= dlr_max; the
published host is the DES-SN DLR assignment matched to DR2): sbw_ne_published (SBW names a host that is not the
published one), dlr_ne_published (DLR names a host that is not the published one), sbw_hostless_published_hosted
(SBW calls the SN hostless although the published host is a DR2 object) and sbw_ne_dlr (both name a host,
different objects). Ellipses as in src/dp2_cases.py (arcsec, angle from +x = West toward North): the SBW one is
the profile's half-light ellipse, the DLR one the ellipse d_DLR is measured on.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import des_lib

des_lib.install()

import dp02_lib as L   # noqa: E402
from dp02_lib import HOSTLESS_CUT, sbw_hostless   # noqa: E402
from dp2_cases import ellipses   # noqa: E402

N_PER_CLASS = 12
SEED = 20261029
ECC = {0: "star", 1: "star?", 2: "galaxy", 3: "galaxy", -9: "no class", -99: "unmatched"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t = pq.read_table(a.associations)
    df = t.to_pandas(ignore_metadata=True)
    attrs = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if not k.startswith(b"pandas")}
    band, profile, dlr_radius = attrs["band"], attrs["profile"], attrs["dlr_radius"]
    s = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    df = df.merge(s[["sn_id"] + [k for k in s.columns if k not in df.columns]], on="sn_id", how="left")
    cut = HOSTLESS_CUT[a.hostless_cut]
    sh, dh = sbw_hostless(df, cut), df.dlr_objectId.to_numpy() < 0
    sb, dl, h = df.sbw_objectId.to_numpy(), df.dlr_objectId.to_numpy(), df.host_objectId.to_numpy()
    classes = {"sbw_ne_published": ~sh & (sb != h), "dlr_ne_published": ~dh & (dl != h),
               "sbw_hostless_published_hosted": sh & (h >= 0), "sbw_ne_dlr": ~sh & ~dh & (sb != dl)}
    rng = np.random.default_rng(SEED)
    picks = []
    for name, sel in classes.items():
        idx = np.flatnonzero(sel)
        idx = rng.choice(idx, size=min(N_PER_CLASS, len(idx)), replace=False)
        picks.append(df.iloc[np.sort(idx)].assign(**{"class": name, "class_n": int(sel.sum())}))
    c = pd.concat(picks, ignore_index=True)

    cand = L.load_candidates(a.candidates, bands=band)
    cand = cand[cand.sn_id.isin(c.sn_id)].drop_duplicates(["sn_id", "objectId"]).reset_index(drop=True)
    cand = pd.concat([cand[["sn_id", "objectId", "coord_ra", "coord_dec"]], ellipses(cand, band, profile, dlr_radius)],
                     axis=1).set_index(["sn_id", "objectId"])
    for name, col in [("sbw", "sbw_objectId"), ("dlr", "dlr_objectId")]:
        pos = cand.reindex(pd.MultiIndex.from_arrays([c.sn_id.to_numpy(), c[col].to_numpy()]))
        c[f"{name}_ra"], c[f"{name}_dec"] = pos.coord_ra.to_numpy(), pos.coord_dec.to_numpy()
        for k in ("a", "b", "phi"):
            c[f"{name}_{k}"] = pos[f"{name}_{k}"].to_numpy()
    c = c.rename(columns={"host_ra": "pub_ra", "host_dec": "pub_dec"})
    c["case_id"] = [f"{k}_{i:02d}" for k, i in zip(c["class"], c.groupby("class").cumcount())]
    c["label"] = [f"{k} | {sid} z {z:.2f} | P(SBW) {p:.2f}, P(pub) {pt:.2f} | flb {f:.2g} | d_DLR {dd:.1f} "
                  f"(pub {pdd:.1f}) | pub {ECC.get(int(e), e)}, {st}"
                  for k, sid, z, p, pt, f, dd, pdd, e, st in zip(
                      c["class"], c.sn_id, c.sn_redshift, c.p_sbw_max, np.nan_to_num(c.p_sbw_true.to_numpy(float)),
                      c.flb_best, c.d_dlr_min, c.host_ddlr_published, c.host_extended_class_coadd, c.true_host_status)]
    cols = (["case_id", "class", "class_n", "sn_id", "sn_ra", "sn_dec", "label", "pub_ra", "pub_dec"]
            + [f"{n}_{k}" for n in ("sbw", "dlr") for k in ("ra", "dec", "a", "b", "phi")]
            + ["p_sbw_max", "p_sbw_true", "p_sbw_dlr", "flb_best", "flb_true", "d_dlr_min", "d_dlr_sbw", "d_dlr_true",
               "host_ddlr_published", "sep_sbw", "sep_dlr", "sep_true", "n_cand_used", "true_host_status",
               "host_extended_class_coadd", "sbw_objectId", "dlr_objectId", "host_objectId", "sn_redshift", "field"])
    c[cols].to_csv(a.output, index=False, float_format="%.10g")


if __name__ == "__main__":
    main()
