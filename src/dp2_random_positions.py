"""Random sky positions for the DP2 hostless control, drawn in each transient's own patch.

For every transient the DP2 Object rows within 20' are read from the HATS partitions (coordinates, tract,
patch only); the transient's patch is that of its nearest object, and its footprint is taken from that
patch's objects (Object rows are primary, so they cover the patch's inner region exactly). ``--per-sn``
positions are drawn uniformly on the sphere inside the objects' RA/Dec bounding box and kept only if they
fall in a 1 arcmin cell holding at least one of the patch's objects, as on DP0.2, so they share the
transient's depth and field (DDF or WFD) and never land on unobserved sky. Several positions per transient
because the clean sample is small (690); each row records its parent transient. Ids are -(100 * k + j + 1)
for the k-th transient (by sn_id) and its j-th position, so they never collide with a DiaObject id.
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

import dp2_lib as D

CELL_ARCMIN = 1.0
READ_RADIUS = 1200.0   # arcsec; > the half-diagonal of a DP2 patch, so the whole patch is read from anywhere in it
_G = {}


def _draw(k):
    sn = _G["sn"]
    sid, ra0, dec0 = int(sn.sn_id.iloc[k]), float(sn.sn_ra.iloc[k]), float(sn.sn_dec.iloc[k])
    paths = _G["boxes"].path.to_numpy()[D.partitions_near(_G["boxes"], ra0, dec0, READ_RADIUS)]
    t = pd.concat([pq.read_table(p, columns=["coord_ra", "coord_dec", "tract", "patch"]).to_pandas() for p in paths],
                  ignore_index=True)
    if t.empty:
        return k, None
    near = np.argmin(np.linalg.norm(D.unit_vectors(t.coord_ra, t.coord_dec) - D.unit_vectors(ra0, dec0), axis=1))
    tract, patch = int(t.tract.iloc[near]), int(t.patch.iloc[near])
    o = t[(t.tract == tract) & (t.patch == patch)]
    ra, dec = o.coord_ra.to_numpy(), o.coord_dec.to_numpy()
    ra = np.where(ra - ra0 > 180, ra - 360, np.where(ra - ra0 < -180, ra + 360, ra))   # unwrap about the SN
    cosd = np.cos(np.deg2rad(np.median(dec)))
    cell = CELL_ARCMIN / 60.0
    r0, r1, d0, d1 = ra.min(), ra.max(), dec.min(), dec.max()
    occ = np.zeros((int((r1 - r0) * cosd / cell) + 1, int((d1 - d0) / cell) + 1), bool)
    occ[((ra - r0) * cosd / cell).astype(int), ((dec - d0) / cell).astype(int)] = True
    rng = np.random.default_rng([_G["seed"], k])
    s0, s1 = np.sin(np.deg2rad(d0)), np.sin(np.deg2rad(d1))
    n, out_r, out_d, have = _G["per_sn"], [], [], 0
    while have < n:
        r = rng.uniform(r0, r1, 4 * n)
        d = np.rad2deg(np.arcsin(rng.uniform(s0, s1, 4 * n)))
        ok = occ[((r - r0) * cosd / cell).astype(int), ((d - d0) / cell).astype(int)]
        out_r.append(r[ok]); out_d.append(d[ok]); have += int(ok.sum())
    r, d = np.concatenate(out_r)[:n] % 360, np.concatenate(out_d)[:n]
    return k, pd.DataFrame({"sn_id": -(100 * k + np.arange(n) + 1), "sn_ra": r, "sn_dec": d,
                            "parent_sn_id": sid, "tract": tract, "patch": patch, "patch_fill": float(occ.mean())})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, help="dp2_sample")
    ap.add_argument("--objects", required=True, help="dp2_object index csv of the HATS partition files")
    ap.add_argument("--per-sn", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261028)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("DP02_WORKERS", "16")))
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    if a.per_sn > 99:
        raise SystemExit("--per-sn must be < 100 (id scheme)")
    t0 = time.time()
    sn = pd.read_parquet(a.sample, columns=["sn_id", "sn_ra", "sn_dec"]).sort_values("sn_id").reset_index(drop=True)
    _G.update(sn=sn, boxes=D.partition_boxes(D.hats_root(a.objects)), per_sn=a.per_sn, seed=a.seed)
    parts = []
    with Pool(a.workers) as pool:
        for k, res in pool.imap_unordered(_draw, range(len(sn))):
            if res is None:
                print(f"no objects within {READ_RADIUS:g}\" of transient {sn.sn_id.iloc[k]}", file=sys.stderr)
            else:
                parts.append(res)
    out = pd.concat(parts, ignore_index=True).sort_values("sn_id").reset_index(drop=True)
    out.to_parquet(a.output, index=False, compression="zstd")
    print(f"wrote {len(out)} random positions for {out.parent_sn_id.nunique()}/{len(sn)} transients in "
          f"{out.groupby(['tract', 'patch']).ngroups} patches (median patch fill {out.patch_fill.median():.2f}) "
          f"in {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
