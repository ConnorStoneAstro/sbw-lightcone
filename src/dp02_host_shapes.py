"""Truth shapes of the DC2 SN host galaxies from cosmoDC2 v1.1.4 (image): bulge and disk
half-light semi-axes, position angle, ellipticities, B/T, stellar masses and magnitudes,
one row per host galaxy_id (the merged-truth SN host_galaxy equals the cosmoDC2 galaxy_id).
Reads only galaxy_id from every row group and the other columns only where a host is."""
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

COLS = (["galaxy_id", "ra", "dec", "redshift", "size_true", "size_minor_true",
         "size_bulge_true", "size_minor_bulge_true", "size_disk_true", "size_minor_disk_true",
         "position_angle_true", "bulge_to_total_ratio_i", "sersic_bulge", "sersic_disk",
         "ellipticity_true", "ellipticity_1_true", "ellipticity_2_true",
         "ellipticity_bulge_true", "ellipticity_1_bulge_true", "ellipticity_2_bulge_true",
         "ellipticity_disk_true", "ellipticity_1_disk_true", "ellipticity_2_disk_true",
         "stellar_mass", "stellar_mass_bulge", "stellar_mass_disk"]
        + [f"mag_true_{b}_lsst" for b in "ugrizy"])

_IDS = None


def _init(ids):
    global _IDS
    _IDS = pa.array(ids, pa.int64())


def _sn_hosts(f):
    t = pq.read_table(f, columns=["host_galaxy", "truth_type"], filters=[("truth_type", "==", 3)])
    return np.unique(t["host_galaxy"].to_numpy())


def _rowgroup(args):
    f, rg = args
    pf = pq.ParquetFile(f)
    gid = pf.read_row_group(rg, columns=["galaxy_id"])["galaxy_id"]
    mask = pc.is_in(gid, value_set=_IDS)
    if not pc.any(mask).as_py():
        return None
    t = pf.read_row_group(rg, columns=COLS).filter(mask)
    return t.to_pandas(ignore_metadata=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, help="DC2 merged truth directory (SN host ids)")
    ap.add_argument("--cosmodc2", required=True, help="csv index of the cosmoDC2 healpix parquet files")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "16")))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    tfiles = sorted(glob.glob(os.path.join(a.truth, "truth_tract*.parquet")))
    with Pool(a.workers) as pool:
        ids = np.unique(np.concatenate(pool.map(_sn_hosts, tfiles)))
    ids = ids[ids > 0]
    print(f"{len(ids)} SN host ids, {time.time() - t0:.0f}s", file=sys.stderr)
    files = pd.read_csv(a.cosmodc2)["path"].tolist()
    tasks = [(f, rg) for f in files for rg in range(pq.ParquetFile(f).num_row_groups)]
    with Pool(a.workers, initializer=_init, initargs=(ids,)) as pool:
        parts = [p for p in pool.imap_unordered(_rowgroup, tasks) if p is not None]
    df = pd.concat(parts, ignore_index=True).drop_duplicates("galaxy_id").sort_values("galaxy_id")
    print(f"{len(df)} of {len(ids)} hosts found in {len(files)} healpix files, {time.time() - t0:.0f}s", file=sys.stderr)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), a.output, compression="zstd")


if __name__ == "__main__":
    main()
