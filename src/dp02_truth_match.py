"""Compare the two DC2 truth matchers (Rubin MatchesTruth vs nearest object within 1 arcsec of the host
truth position): how often they name a different host object, and SBW / DLR recovery where they agree."""
from __future__ import annotations

import argparse
import json

import numpy as np
import pyarrow.parquet as pq

from dp02_lib import HOSTLESS_CUT, matrix_cells, wilson


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--hostless-cut", required=True, choices=sorted(HOSTLESS_CUT))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    cut = HOSTLESS_CUT[a.hostless_cut]
    precut = json.loads(pq.read_schema(a.sample).metadata[b"truth_match_precut"])
    df = pq.read_table(a.associations).to_pandas(ignore_metadata=True)
    mt, nn = df.host_objectId_mt.to_numpy(), df.host_objectId_nn.to_numpy()
    both = (mt >= 0) & (nn >= 0)
    same = both & (mt == nn)
    res = {"before_matcher_cut": precut,
           "sample": {"n": len(df), "both": int(both.sum()), "mt_only": int(((mt >= 0) & (nn < 0)).sum()),
                      "nn_only": int(((mt < 0) & (nn >= 0)).sum()), "both_differ": int((both & ~same).sum()),
                      "median_match_sep_arcsec_nn": float(np.nanmedian(df.match_sep_arcsec))}}
    for name, sel in [("matchers_agree", same), ("matchers_differ", both & ~same)]:
        d = df[sel]
        r = {"n": len(d)}
        for m in ("sbw", "dlr"):
            for tn, tc in [("mt", "host_objectId_mt"), ("nn", "host_objectId_nn")]:
                k = int(((d[f"{m}_objectId"] == d[tc]) & (d[tc] >= 0)).sum())
                p, lo, hi = wilson(k, len(d))
                r[f"{m}_recovery_vs_{tn}"] = float(p); r[f"{m}_recovery_vs_{tn}_ci68"] = [float(lo), float(hi)]
        r["agreement"] = float(d.agree.mean())
        res[name] = r
    res["matrix_own_truth"] = matrix_cells(df, cut)
    res["matrix_common_vs_mt"] = matrix_cells(df[both], cut, "host_objectId_mt")
    res["matrix_common_vs_nn"] = matrix_cells(df[both], cut, "host_objectId_nn")
    res["sbw_hostless_cut"] = cut
    with open(a.output, "w") as f:
        json.dump(res, f, indent=2)


if __name__ == "__main__":
    main()
