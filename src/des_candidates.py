"""Every DES DR2 object within the search radius of each position (des_sample SNe or des_random_positions), one
row per (position, object) pair in the candidate schema of src/des_lib.py, from the 2-arcmin des_dr2_objects
cache. Refuses to run if a position's cone is not inside the cached cone of its SN (parent_sn_id for random
positions)."""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

import des_lib as D
from dp02_lib import RADIUS


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--positions", required=True, help="table with sn_id, sn_ra, sn_dec (and parent_sn_id)")
    ap.add_argument("--sample", required=True, help="des_sample (the SN positions the cache was fetched around)")
    ap.add_argument("--objects", required=True)
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    r = RADIUS[a.search_radius]
    t = pq.read_table(a.objects)
    cache_r = float(t.schema.metadata[b"radius_arcsec"])
    obj = D.candidate_columns(t.to_pandas(ignore_metadata=True))
    pos = pq.read_table(a.positions).to_pandas(ignore_metadata=True)
    sn = pq.read_table(a.sample, columns=["sn_id", "sn_ra", "sn_dec"]).to_pandas().set_index("sn_id")
    parent = pos.parent_sn_id if "parent_sn_id" in pos else pos.sn_id
    pr = sn.reindex(parent.to_numpy())
    off = D.arcsec_from_chord(np.linalg.norm(D.unit_vectors(pos.sn_ra, pos.sn_dec)
                                             - D.unit_vectors(pr.sn_ra, pr.sn_dec), axis=1))
    if not (off + r <= cache_r + 1e-6).all():
        raise SystemExit(f"{int((off + r > cache_r).sum())} cones leave the {cache_r:g}\" cache")
    ou = D.unit_vectors(obj.coord_ra, obj.coord_dec)
    tree = cKDTree(ou)
    pv = D.unit_vectors(pos.sn_ra, pos.sn_dec)
    hits = tree.query_ball_point(pv, D.chord(r))
    n = np.array([len(k) for k in hits])
    rows = np.repeat(np.arange(len(pos)), n)
    idx = np.concatenate([np.asarray(k, int) for k in hits]) if n.sum() else np.zeros(0, int)
    head = pd.DataFrame({"sn_id": pos.sn_id.to_numpy(np.int64)[rows],
                         "sep_arcsec": D.arcsec_from_chord(np.linalg.norm(ou[idx] - pv[rows], axis=1))})
    out = pd.concat([head, obj.iloc[idx].reset_index(drop=True)], axis=1)
    out = out.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)
    tab = pa.Table.from_pandas(out, preserve_index=False)
    tab = tab.replace_schema_metadata({b"search_radius": a.search_radius.encode(), b"n_positions": str(len(pos)).encode(),
                                       b"n_without_candidates": str(int((n == 0).sum())).encode()})
    pq.write_table(tab, a.output, compression="zstd")
    print(f"{len(out)} pairs for {len(pos)} positions within {r:g}\" (median {np.median(n):.0f}, none for "
          f"{int((n == 0).sum())})", file=sys.stderr)


if __name__ == "__main__":
    main()
