"""Extract every DP0.2 primary object within the search radius of any DC2 SN.

Reads the SN positions from the DESC merged truth (truth_type == 3) and walks
the 157 DP0.2 objectTable tract files in parallel, keeping one row per
(SN, object) pair with the columns the SBW and DLR models need.  Output is a
single parquet file; everything downstream reads it instead of the 1 TB
Object table.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy.spatial import cKDTree

BANDS = "ugrizy"
PER_BAND = ["cModelFlux", "bdReB", "bdReD", "bdFluxB", "bdFluxD",
            "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF", "kronRad", "extendedness"]
GENERAL = ["objectId", "coord_ra", "coord_dec", "detect_isPrimary", "refExtendedness", "refBand",
           "tract", "patch", "parentObjectId"]
RADIUS = {"r15": 15.0, "r30": 30.0, "r60": 60.0}


def unit_vectors(ra_deg, dec_deg):
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def load_sne(truth_dir):
    frames = []
    for f in sorted(glob.glob(os.path.join(truth_dir, "truth_tract*.parquet"))):
        t = pq.read_table(f, columns=["id", "truth_type", "ra", "dec"],
                          filters=[("truth_type", "==", 3)])
        frames.append(t.select(["id", "ra", "dec"]).to_pandas())
    sn = pd.concat(frames, ignore_index=True).drop_duplicates("id")
    return sn


def one_tract(args):
    path, sn_ra, sn_dec, sn_id, radius_arcsec = args
    cols = GENERAL + [f"{b}_{c}" for b in BANDS for c in PER_BAND]
    pf = pq.ParquetFile(path)
    have = set(pf.schema_arrow.names)
    cols = [c for c in cols if c in have]
    t = pq.read_table(path, columns=cols, filters=[("detect_isPrimary", "==", True)])
    if t.num_rows == 0:
        return None
    ora = t["coord_ra"].to_numpy(zero_copy_only=False)
    odec = t["coord_dec"].to_numpy(zero_copy_only=False)
    pad = radius_arcsec / 3600.0 + 0.01
    # SNe that can have a candidate in this tract: inside the object bounding box plus padding
    cosd = np.cos(np.deg2rad(np.median(odec)))
    box = ((sn_dec >= odec.min() - pad) & (sn_dec <= odec.max() + pad)
           & (sn_ra >= ora.min() - pad / cosd) & (sn_ra <= ora.max() + pad / cosd))
    if not box.any():
        return None
    sra, sdec, sid = sn_ra[box], sn_dec[box], sn_id[box]
    tree = cKDTree(unit_vectors(ora, odec))
    chord = 2.0 * np.sin(np.deg2rad(radius_arcsec / 3600.0) / 2.0)
    hits = tree.query_ball_point(unit_vectors(sra, sdec), r=chord)
    n = np.fromiter((len(h) for h in hits), int, len(hits))
    if n.sum() == 0:
        return None
    obj_idx = np.concatenate([np.asarray(h, int) for h in hits if len(h)])
    sn_rep = np.repeat(np.arange(len(hits)), n)
    sep = np.rad2deg(2 * np.arcsin(np.linalg.norm(
        unit_vectors(sra[sn_rep], sdec[sn_rep]) - unit_vectors(ora[obj_idx], odec[obj_idx]), axis=1) / 2)) * 3600
    sub = t.take(pa.array(obj_idx)).drop_columns(["detect_isPrimary"]).to_pandas()
    out = pd.DataFrame({"sn_id": sid[sn_rep].astype(np.int64), "sep_arcsec": sep.astype(np.float32)})
    for c in sub.columns:
        v = sub[c]
        if c in ("coord_ra", "coord_dec"):
            out[c] = v.to_numpy(dtype=np.float64)
        elif v.dtype.kind == "f":
            out[c] = v.to_numpy(dtype=np.float32)
        else:
            out[c] = v.to_numpy()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, help="directory of DC2 merged truth tract parquet files")
    ap.add_argument("--objects", required=True, help="directory of DP0.2 objectTable tract parquet files")
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "16")))
    ap.add_argument("--max-tracts", type=int, default=None, help="probe: only the first N tract files")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)

    t0 = time.time()
    sn = load_sne(a.truth)
    print(f"{len(sn)} SNe from truth", file=sys.stderr)
    files = sorted(glob.glob(os.path.join(a.objects, "objectTable_tract_*.parq")))
    if a.max_tracts:
        files = files[: a.max_tracts]
    r = RADIUS[a.search_radius]
    jobs = [(f, sn["ra"].to_numpy(), sn["dec"].to_numpy(), sn["id"].to_numpy(), r) for f in files]
    parts = []
    with Pool(a.workers) as pool:
        for i, part in enumerate(pool.imap_unordered(one_tract, jobs)):
            if part is not None:
                parts.append(part)
            if (i + 1) % 10 == 0:
                print(f"{i + 1}/{len(files)} tracts, {sum(len(p) for p in parts)} pairs, {time.time() - t0:.0f}s",
                      file=sys.stderr)
    cand = pd.concat(parts, ignore_index=True)
    cand = cand.drop_duplicates(["sn_id", "objectId"]).sort_values(["sn_id", "sep_arcsec"]).reset_index(drop=True)
    cand.attrs = {}
    table = pa.Table.from_pandas(cand, preserve_index=False)
    table = table.replace_schema_metadata({b"search_radius_arcsec": str(r).encode(),
                                           b"n_sne_with_candidates": str(cand["sn_id"].nunique()).encode()})
    pq.write_table(table, a.output, compression="zstd")
    print(f"wrote {len(cand)} pairs for {cand['sn_id'].nunique()} SNe to {a.output} in {time.time() - t0:.0f}s",
          file=sys.stderr)


if __name__ == "__main__":
    main()
