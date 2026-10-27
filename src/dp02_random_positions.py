"""One random sky position per DC2 truth SN, drawn uniformly in that SN's tract.

The tract's footprint is taken from its DP0.2 objectTable file: the primary
objects (detect_isPrimary) cover exactly the tract's inner region, so a
position is drawn uniformly on the sphere inside their RA/Dec bounding box
(RA uniform, sin(dec) uniform) and kept only if it falls in a 1 arcmin cell
that holds at least one primary object; the rest are redrawn.  This keeps
random positions off the unobserved corners of edge tracts, where they would
be trivially hostless.  SNe in tracts without an objectTable file get no
random position.  The random ids are the negated truth SN ids, so they never
collide with a real sn_id.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

CELL_ARCMIN = 1.0


def _sne(f):
    t = pq.read_table(f, columns=["id", "truth_type", "tract"], filters=[("truth_type", "==", 3)])
    return t.select(["id", "tract"]).to_pandas(ignore_metadata=True)


def _draw(args):
    path, tract, n, seed = args
    t = pq.read_table(path, columns=["coord_ra", "coord_dec", "detect_isPrimary"],
                      filters=[("detect_isPrimary", "==", True)])
    ra = t["coord_ra"].to_numpy(zero_copy_only=False)
    dec = t["coord_dec"].to_numpy(zero_copy_only=False)
    if len(ra) == 0:
        return tract, None
    cosd = np.cos(np.deg2rad(np.median(dec)))
    cell = CELL_ARCMIN / 60.0
    ra0, ra1, d0, d1 = ra.min(), ra.max(), dec.min(), dec.max()
    nx = int(np.ceil((ra1 - ra0) * cosd / cell)) + 1
    ny = int(np.ceil((d1 - d0) / cell)) + 1
    occ = np.zeros((nx, ny), bool)
    occ[((ra - ra0) * cosd / cell).astype(int), ((dec - d0) / cell).astype(int)] = True
    rng = np.random.default_rng([seed, int(tract)])
    s0, s1 = np.sin(np.deg2rad(d0)), np.sin(np.deg2rad(d1))
    out_ra, out_dec, have = [], [], 0
    while have < n:
        m = 2 * (n - have) + 100
        r = rng.uniform(ra0, ra1, m)
        d = np.rad2deg(np.arcsin(rng.uniform(s0, s1, m)))
        ok = occ[((r - ra0) * cosd / cell).astype(int), ((d - d0) / cell).astype(int)]
        out_ra.append(r[ok]); out_dec.append(d[ok]); have += int(ok.sum())
    return tract, (np.concatenate(out_ra)[:n], np.concatenate(out_dec)[:n], float(occ.mean()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, help="directory of DC2 merged truth tract parquet files")
    ap.add_argument("--objects", required=True, help="csv index of the DP0.2 objectTable tract files")
    ap.add_argument("--seed", type=int, default=20261027)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "16")))
    ap.add_argument("--max-tracts", type=int, default=None, help="probe: only the first N tracts with SNe")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    tfiles = sorted(glob.glob(os.path.join(a.truth, "truth_tract*.parquet")))
    with Pool(a.workers) as pool:
        sn = pd.concat(pool.map(_sne, tfiles), ignore_index=True).drop_duplicates("id")
    sn = sn.sort_values("id").reset_index(drop=True)
    files = {int(re.search(r"objectTable_tract_(\d+)_", p).group(1)): p
             for p in pd.read_csv(a.objects)["path"]}
    counts = sn.groupby("tract").size()
    tracts = [int(t) for t in counts.index if int(t) in files]
    if a.max_tracts:
        tracts = tracts[: a.max_tracts]
    print(f"{len(sn)} SNe in {len(counts)} tracts, {len(tracts)} with an objectTable file", file=sys.stderr)
    jobs = [(files[t], t, int(counts.loc[t]), a.seed) for t in tracts]
    parts = []
    with Pool(a.workers) as pool:
        for tract, res in pool.imap_unordered(_draw, jobs):
            if res is None:
                continue
            r, d, filled = res
            ids = sn.loc[sn.tract == tract, "id"].to_numpy(np.int64)
            parts.append(pd.DataFrame({"sn_id": -ids, "sn_ra": r, "sn_dec": d,
                                       "tract": np.full(len(ids), tract, np.int32)}))
            if filled < 0.5:
                print(f"tract {tract}: only {filled:.0%} of the bounding box is covered", file=sys.stderr)
    out = pd.concat(parts, ignore_index=True).sort_values("sn_id").reset_index(drop=True)
    assert (out.sn_id < 0).all()
    out.to_parquet(a.output, index=False, compression="zstd")
    print(f"wrote {len(out)} random positions in {out.tract.nunique()} tracts "
          f"({len(sn) - len(out)} SNe without one) in {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
