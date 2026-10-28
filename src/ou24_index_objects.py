"""Inventory the OU2024 w_2024_22 objectTable_tract files (step3b_wfd_* runs).

One row per file: run, tract, path, size, rows, primary rows, the patches it
holds, and the RA/Dec bounding box of its primary objects. Also writes a
per-tract summary and a column check against the DP0.2 needs. The output csv
(``path`` column) is the objects index read by ``ou24_candidates.py``.

Several step3b runs cover the same tract. ``--rule`` decides which files enter
the index: ``union`` keeps every file (correct when the runs split a tract by
patch), ``latest`` keeps the highest-numbered run per tract.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = "/global/cfs/cdirs/lsst/production/gen3/roman-desc-sims/repo/u/descdm"
DP02_PER_BAND = ["cModelFlux", "bdReB", "bdReD", "bdFluxB", "bdFluxD",
                 "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF", "kronRad", "extendedness"]
DP02_GENERAL = ["objectId", "coord_ra", "coord_dec", "detect_isPrimary", "refExtendedness", "refBand",
                "tract", "patch", "parentObjectId"]
EXTRA_PER_BAND = ["bdE1", "bdE2", "bdChi2", "gaapOptimalFlux", "gaap1p0Flux", "calibFlux", "psfFlux",
                  "sizeExtendedness", "blendedness", "deblend_blendedness", "hsmShapeRegauss_e1",
                  "ixxDebiasedPSF", "kronFlux", "free_cModelFlux", "cModel_flag", "i_flag", "fwhm"]
EXTRA_GENERAL = ["refSizeExtendedness", "deblend_nChild", "deblend_failed", "deblend_skipped",
                 "detect_isDeblendedModelSource", "detect_isIsolated", "detect_fromBlend", "ebv", "refFwhm",
                 "sky_object", "footprintArea"]


def one_file(path):
    m = re.search(r"(step3b_wfd_(\d+)_w_2024_22)/.*objectTable_tract/(\d+)/", path)
    pf = pq.ParquetFile(path)
    t = pq.read_table(path, columns=["patch", "detect_isPrimary", "coord_ra", "coord_dec"])
    prim = t["detect_isPrimary"].to_numpy(zero_copy_only=False).astype(bool)
    patch = t["patch"].to_numpy(zero_copy_only=False)
    ra = t["coord_ra"].to_numpy(zero_copy_only=False)[prim]
    dec = t["coord_dec"].to_numpy(zero_copy_only=False)[prim]
    return dict(run=m.group(1), run_number=int(m.group(2)), tract=int(m.group(3)), path=path,
                size_bytes=os.path.getsize(path), mtime=int(os.path.getmtime(path)),
                n_rows=pf.metadata.num_rows, n_primary=int(prim.sum()),
                patches=" ".join(map(str, sorted(set(patch.tolist())))),
                ra_min=float(ra.min()), ra_max=float(ra.max()), dec_min=float(dec.min()), dec_max=float(dec.max()))


def bde1_check(path, band="r", nmax=200000):
    """Compare bdE1/bdE2 with the ellipticity implied by the CModel and the
    PSF-subtracted HSM moments. The v23 functor bug gave bdE1 = xx - yy/(xx+yy)."""
    cols = [f"{band}_{c}" for c in ("bdE1", "bdE2", "bdReD", "bdReB", "bdFluxD", "bdFluxB",
                                     "ixx", "iyy", "ixy", "ixxPSF", "iyyPSF", "ixyPSF")] + ["detect_isPrimary", "refExtendedness"]
    d = pq.read_table(path, columns=cols).to_pandas(ignore_metadata=True)
    d = d[d.detect_isPrimary & (d.refExtendedness > 0.5)].head(nmax)
    g = lambda k: d[f"{band}_{k}"].to_numpy(float)
    mxx, myy, mxy = g("ixx") - g("ixxPSF"), g("iyy") - g("iyyPSF"), g("ixy") - g("ixyPSF")
    e1m = (mxx - myy) / (mxx + myy)
    e2m = 2 * mxy / (mxx + myy)
    e1, e2 = g("bdE1"), g("bdE2")
    ok = np.isfinite(e1) & np.isfinite(e1m) & np.isfinite(e2) & (mxx > 0) & (myy > 0) & (np.abs(e1m) < 1) \
        & (g("bdReD") > 2)                                       # resolved, so the shapes should agree
    e1, e2, e1m, e2m = e1[ok], e2[ok], e1m[ok], e2m[ok]
    # distortion (a^2-b^2)/(a^2+b^2) vs shear-like (a-b)/(a+b): compare sign agreement and correlation
    out = dict(n=int(ok.sum()),
               e1_range=[float(np.nanpercentile(e1, 1)), float(np.nanpercentile(e1, 99))],
               e2_range=[float(np.nanpercentile(e2, 1)), float(np.nanpercentile(e2, 99))],
               frac_abs_e1_gt1=float(np.mean(np.abs(e1) > 1)),
               corr_e1=float(np.corrcoef(e1, e1m)[0, 1]), corr_e2=float(np.corrcoef(e2, e2m)[0, 1]),
               sign_agree_e1=float(np.mean(np.sign(e1) == np.sign(e1m))),
               sign_agree_e2=float(np.mean(np.sign(e2) == np.sign(e2m))),
               median_ratio_e1=float(np.median(e1 / e1m)),
               median_ratio_e2=float(np.median(e2 / e2m)))
    # the v23 bug signature: bdE1 tracks xx - yy/(xx+yy) of the model ellipse, i.e. grows with size
    out["corr_abs_e1_with_reD"] = float(np.corrcoef(np.abs(e1), g("bdReD")[ok])[0, 1])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--rule", required=True, choices=["union", "latest"])
    ap.add_argument("--band", default="r")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output", required=True, help="index csv (one row per selected file)")
    ap.add_argument("--all-output", help="csv of every file found, selected or not")
    ap.add_argument("--report", help="json with per-tract summary, column check and the bdE1 test")
    a = ap.parse_args(argv)

    files = sorted(glob.glob(os.path.join(a.repo, "step3b_wfd_*_w_2024_22", "*", "objectTable_tract", "*", "*.parq")))
    print(f"{len(files)} objectTable files", file=sys.stderr)
    with Pool(a.workers) as pool:
        rows = pool.map(one_file, files)
    df = pd.DataFrame(rows).sort_values(["tract", "run_number"]).reset_index(drop=True)
    if a.all_output:
        df.to_csv(a.all_output, index=False)
    # do the runs of one tract share patches?
    overlap = {}
    for tr, g in df.groupby("tract"):
        sets = [set(p.split()) for p in g.patches]
        overlap[int(tr)] = int(sum(len(s1 & s2) for i, s1 in enumerate(sets) for s2 in sets[i + 1:]))
    sel = df if a.rule == "union" else df.sort_values("run_number").groupby("tract").tail(1)
    sel = sel.sort_values(["tract", "run_number"]).reset_index(drop=True)
    sel.to_csv(a.output, index=False)
    print(f"selected {len(sel)} files ({a.rule}) for {sel.tract.nunique()} tracts", file=sys.stderr)

    if a.report:
        names = set(pq.ParquetFile(files[0]).schema_arrow.names)
        need = DP02_GENERAL + [f"{b}_{c}" for b in "ugrizy" for c in DP02_PER_BAND]
        extra = EXTRA_GENERAL + [f"r_{c}" for c in EXTRA_PER_BAND]
        tr = sel.groupby("tract").agg(n_files=("path", "size"), runs=("run", " ".join), n_rows=("n_rows", "sum"),
                                      n_primary=("n_primary", "sum"), ra_min=("ra_min", "min"), ra_max=("ra_max", "max"),
                                      dec_min=("dec_min", "min"), dec_max=("dec_max", "max"))
        tr["patch_overlap_between_runs"] = [overlap[t] for t in tr.index]
        dec0 = np.deg2rad((tr.dec_min + tr.dec_max) / 2)
        tr["area_bbox_deg2"] = (tr.ra_max - tr.ra_min) * np.cos(dec0) * (tr.dec_max - tr.dec_min)
        checks = {}
        for f in sel.groupby("tract").head(1).path.tolist()[:3]:
            checks[os.path.basename(f)] = bde1_check(f, a.band)
        rep = dict(n_files_found=len(files), n_files_selected=len(sel), rule=a.rule,
                   n_tracts=int(sel.tract.nunique()), n_columns=len(names),
                   missing_dp02_columns=[c for c in need if c not in names],
                   extra_columns_present={c: c in names for c in extra},
                   total_rows=int(sel.n_rows.sum()), total_primary=int(sel.n_primary.sum()),
                   total_size_gb=round(sel.size_bytes.sum() / 1e9, 1),
                   total_bbox_area_deg2=round(float(tr.area_bbox_deg2.sum()), 2),
                   per_tract=json.loads(tr.reset_index().to_json(orient="records")),
                   bde1_check=checks)
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)


if __name__ == "__main__":
    main()
