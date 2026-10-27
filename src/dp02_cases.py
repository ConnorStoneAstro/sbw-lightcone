"""Random disagreement cases for the cutout gallery: for each class, a seeded draw of SNe with the
positions of every marker (truth host position, both truth matchers' objects, the SBW most probable host and
the DLR host)."""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from dp02_lib import P_TRUE, p_true

N_PER_CLASS = 12
SEED = 20261027


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--p-threshold", required=True, choices=sorted(P_TRUE))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    df = pq.read_table(a.associations).to_pandas(ignore_metadata=True)
    t, mt, nn = df.host_objectId.to_numpy(), df.host_objectId_mt.to_numpy(), df.host_objectId_nn.to_numpy()
    s, d = df.sbw_objectId.to_numpy(), df.dlr_objectId.to_numpy()
    pt, st = p_true(df), df.true_host_status.to_numpy()
    thr = P_TRUE[a.p_threshold]
    classes = {"truth_matchers_disagree": mt != nn,
               "sbw_excludes_true_host": (st == "kept") & (pt < thr),
               "dlr_wrong_or_hostless": d != t,
               "sbw_most_probable_ne_dlr": s != d,
               "true_host_cut_by_object_table": np.isin(st, ["point_source", "unusable", "quality"])}
    rng = np.random.default_rng(SEED)
    picks = []
    for name, sel in classes.items():
        idx = np.flatnonzero(sel)
        idx = rng.choice(idx, size=min(N_PER_CLASS, len(idx)), replace=False)
        picks.append(df.iloc[np.sort(idx)].assign(**{"class": name, "class_n": int(sel.sum())}))
    c = pd.concat(picks, ignore_index=True)

    cand = pq.read_table(a.candidates, columns=["sn_id", "objectId", "coord_ra", "coord_dec"],
                         filters=[("sn_id", "in", c.sn_id.unique().tolist())]).to_pandas(ignore_metadata=True)
    cand = cand.drop_duplicates(["sn_id", "objectId"]).set_index(["sn_id", "objectId"])
    for name, col in [("truth_match", "host_objectId_mt"), ("truth_nearest", "host_objectId_nn"),
                      ("sbw", "sbw_objectId"), ("dlr", "dlr_objectId")]:
        key = pd.MultiIndex.from_arrays([c.sn_id.to_numpy(), c[col].to_numpy()])
        pos = cand.reindex(key)
        c[f"{name}_ra"], c[f"{name}_dec"] = pos.coord_ra.to_numpy(), pos.coord_dec.to_numpy()
    c["truth_pos_ra"], c["truth_pos_dec"] = c.host_ra, c.host_dec
    c["case_id"] = [f"{k}_{i:02d}" for k, i in zip(c["class"], c.groupby("class").cumcount())]
    c["p_sbw_true"] = p_true(c)
    c["label"] = [f"{k} | SN {sid} | P(true) {pt_:.3f}, P(most prob.) {p:.2f} | {stt} | d_DLR {dd:.1f} | sep {sp:.1f}\""
                  for k, sid, pt_, p, stt, dd, sp in zip(c["class"], c.sn_id, c.p_sbw_true, c.p_sbw_max,
                                                         c.true_host_status, c.d_dlr_min, c.sep_truth_arcsec)]
    cols = ["case_id", "class", "class_n", "sn_id", "sn_ra", "sn_dec", "label",
            "truth_pos_ra", "truth_pos_dec", "truth_match_ra", "truth_match_dec",
            "truth_nearest_ra", "truth_nearest_dec", "sbw_ra", "sbw_dec", "dlr_ra", "dlr_dec",
            "p_sbw_max", "p_sbw_true", "true_host_status", "d_dlr_min", "sep_truth_arcsec",
            "host_objectId_mt", "host_objectId_nn", "sbw_objectId", "dlr_objectId"]
    c[cols].to_csv(a.output, index=False, float_format="%.7f")


if __name__ == "__main__":
    main()
