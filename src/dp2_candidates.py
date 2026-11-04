"""Every DP2 Object within the search radius of each transient (or random position), one row per pair.

Objects come from the NERSC HATS copy of the DP2 Object table (object_collection/object_lc, 786M rows in
8840 healpix partitions; Object rows are all primary). Instead of an LSDB crossmatch, every partition whose
coord bounding box (per_partition_statistics.parquet), padded by the radius, contains a position is read
whole, so there is no margin cache to get wrong: a position near a partition edge simply reads both
partitions, which are disjoint. Each partition is read twice, first coord_ra / coord_dec to match, then the
candidate columns of the matched rows only. Columns are those of src/dp2_lib.py (both MultiProFit fits, per
band model / cModel fluxes, HSM moments, extendedness, blendedness, cModel_flag, deblend and parent columns;
deblend_parentNChild is written as parent_nChild like on DP0.2).

"""
from __future__ import annotations

import argparse
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

import dp2_lib as D
from dp02_lib import RADIUS

_G = {}


def _match(job):
    path, idx = job
    sn = _G["sn"]
    t = pq.read_table(path, columns=["coord_ra", "coord_dec"])
    ra, dec = t["coord_ra"].to_numpy(), t["coord_dec"].to_numpy()
    if len(ra) == 0:
        return None
    tree = cKDTree(D.unit_vectors(ra, dec))
    hits = tree.query_ball_point(D.unit_vectors(sn.sn_ra.to_numpy()[idx], sn.sn_dec.to_numpy()[idx]),
                                 D.chord(_G["radius"]))
    pairs = [(i, j) for i, js in zip(idx, hits) for j in js]
    if not pairs:
        return None
    si, oj = np.array(pairs).T
    rows = np.unique(oj)
    cand = pq.read_table(path, columns=_G["columns"]).take(rows).to_pandas(ignore_metadata=True)
    pos = np.searchsorted(rows, oj)
    out = cand.iloc[pos].reset_index(drop=True)
    out.insert(0, "sn_id", sn.sn_id.to_numpy()[si])
    a, b = D.unit_vectors(out.coord_ra, out.coord_dec), D.unit_vectors(sn.sn_ra.to_numpy()[si], sn.sn_dec.to_numpy()[si])
    out.insert(4, "sep_arcsec", np.degrees(2 * np.arcsin(np.clip(np.linalg.norm(a - b, axis=1) / 2, 0, 1))) * 3600)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--positions", required=True, help="parquet with sn_id, sn_ra, sn_dec (dp2_sample or random positions)")
    ap.add_argument("--catalogue", choices=["dp2_object"], default="dp2_object",
                    help="fixed; the Legacy DR10 alternative was dropped on 4 Nov 2026")
    ap.add_argument("--objects", required=True, help="dp2_object index csv of the HATS partition files")
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "16")))
    ap.add_argument("--record-decisions", action="store_true", help="write the active decisions into the parquet metadata")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    sn = pd.read_parquet(a.positions, columns=["sn_id", "sn_ra", "sn_dec"]).reset_index(drop=True)
    radius = RADIUS[a.search_radius]
    root = D.hats_root(a.objects)
    boxes = D.partition_boxes(root)
    listed = set(pd.read_csv(a.objects)["path"])
    missing = set(boxes.path) - listed
    if missing:
        raise SystemExit(f"{len(missing)} HATS partitions are not in {a.objects}, e.g. {sorted(missing)[0]}")
    jobs = {}
    for i, (r, d) in enumerate(zip(sn.sn_ra, sn.sn_dec)):
        for k in D.partitions_near(boxes, r, d, radius):
            jobs.setdefault(boxes.path.iloc[k], []).append(i)
    print(f"{len(sn)} positions, {len(jobs)} partitions to read", file=sys.stderr)
    _G.update(sn=sn, radius=radius, columns=D.object_columns())
    parts = []
    with Pool(a.workers) as pool:
        for out in pool.imap_unordered(_match, [(p, np.array(ix)) for p, ix in jobs.items()]):
            if out is not None:
                parts.append(out)
    cand = pd.concat(parts, ignore_index=True).rename(columns=D.RENAME)
    write(cand, a.output, len(sn), radius, t0, {"dp2_candidate_catalogue": a.catalogue} if a.record_decisions else None)


def write(cand, output, n_pos, radius, t0, decisions=None):
    cand = cand.sort_values(["sn_id", "sep_arcsec"], kind="stable").reset_index(drop=True)
    for c in ("sn_id", "objectId", "parentObjectId"):
        cand[c] = cand[c].fillna(0).astype(np.int64)
    if decisions:
        import json as _json
        import pyarrow as _pa
        import pyarrow.parquet as _pq
        tab = _pa.Table.from_pandas(cand, preserve_index=False)
        tab = tab.replace_schema_metadata({**(tab.schema.metadata or {}),
                                           b"decisions": _json.dumps(decisions).encode()})
        _pq.write_table(tab, output, compression="zstd")
    else:
        cand.to_parquet(output, index=False, compression="zstd")
    n_with = cand.sn_id.nunique()
    print(f"wrote {len(cand)} pairs for {n_with}/{n_pos} positions "
          f"(median {cand.groupby('sn_id').size().median():.0f} objects within {radius:g}\") in {time.time() - t0:.0f}s",
          file=sys.stderr)


if __name__ == "__main__":
    main()
