"""OU2024 transient table: every transient in the 22 w_2024_22 tracts with its true host and the host's
skyCatalog truth (position, magnitudes, bulge/disk half-light radii, shapes, bulge/disk/knot light
fractions, stellar mass).

Transient truth (``--truth-index``, files picked by name):
* skyCatalog ``snana_<hpix>.parquet`` (what imSim injected): id, ra, dec, host_id, z_CMB, model_name, mjds.
* SNANA dump ``OU24_TRANSIENT_TRUTH_HEAD.DAT.gz``: CID (= skyCatalog id), GENTYPE, GALID, GALSNSEP,
  GALSNDDLR, LOGMASS_TRUE, a_rot, PEAKMAG_r.
* Per-class ``ROMAN+LSST_LARGE_<class>/*_HEAD.FITS.gz``: SIM_HOSTLIB_GALID (true host) and SNANA's
  HOSTGAL_* match (HOSTGAL_OBJID is the true host whenever non-zero, so it is not an independent DLR host).
The true host is GALID (= SIM_HOSTLIB_GALID). The skyCatalog host_id is 0 for ~2.4% of transients.

Host truth (``--galaxy-index``): skyCatalog ``galaxy_<hpix>.parquet`` (positions, half-light radii,
reduced-shear ellipticities, stellar mass), ``galaxy_flux_<hpix>.parquet`` (LSST fluxes, photons/s/cm^2,
turned into AB mags with the 00README zero points) and ``galaxy_sed_<hpix>.hdf5`` (disk, bulge and knot
SEDs, integrated through the ``total_<b>.dat`` throughputs for the component light fractions; knots share
the disk shape, so B/T is the bulge fraction). galaxy_id = <hpix><1 digit><9 digits>.

Shapes: the half-light radii are those of the galsim Sersic profile before the (area-preserving) shear,
so they are determinant radii sqrt(a b): semi-major = hlr / sqrt(q). q = (1-|g|)/(1+|g|),
beta = 0.5 atan2(g2, g1) in degrees. beta equals SNANA's a_rot and is the major-axis angle from +x (West)
toward +y (North) in the sbw_host.core.sky_to_plane frame: SN offsets of elongated hosts align with it
(<cos 2(psi - beta)> = 0.51 for q < 0.5) and so do the pipeline moments of the matched object
(see ou24_pa_check).

A transient is kept when a primary object lies within ``--coverage-radius`` arcsec; its ``tract`` is
that object's tract.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from astropy.io import fits
from scipy.spatial import cKDTree

SNANA_COLS = ["id", "ra", "dec", "host_id", "host_ra", "host_dec", "host_sn_sep", "z_CMB", "gentype", "model_name",
              "start_mjd", "end_mjd", "peak_mjd", "peak_mag_g", "peak_mag_i"]
DAT_COLS = ["CID", "GENTYPE", "GALID", "GALSNSEP", "GALSNDDLR", "LOGMASS_TRUE", "LOG_SFR", "a_rot", "PEAKMAG_r"]
HEAD_COLS = ["SNID", "SIM_HOSTLIB_GALID", "HOSTGAL_OBJID", "HOSTGAL_DDLR", "HOSTGAL_SNSEP", "HOSTGAL_NMATCH",
             "HOSTGAL_MAG_r", "HOSTGAL_LOGMASS", "SIM_HOSTLIB(LOGMASS_TRUE)"]
ZP = {"u": 13.154216481475633, "g": 14.562127424949786, "r": 14.354587987629897,
      "i": 14.064951243337182, "z": 13.672316495190309, "y": 12.881577415640477}
GAL_COLS = ["galaxy_id", "ra", "dec", "redshift", "spheroidHalfLightRadiusArcsec", "diskHalfLightRadiusArcsec",
            "diskEllipticity1", "diskEllipticity2", "spheroidEllipticity1", "spheroidEllipticity2",
            "um_source_galaxy_obs_sm"]
RENAME = {"galaxy_id": "host_id", "ra": "host_ra", "dec": "host_dec", "redshift": "host_redshift",
          "spheroidHalfLightRadiusArcsec": "host_hlr_bulge", "diskHalfLightRadiusArcsec": "host_hlr_disk",
          "diskEllipticity1": "host_g1_disk", "diskEllipticity2": "host_g2_disk",
          "spheroidEllipticity1": "host_g1_bulge", "spheroidEllipticity2": "host_g2_bulge",
          "um_source_galaxy_obs_sm": "host_stellar_mass"}


def unit_vectors(ra_deg, dec_deg):
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _objects(path):
    t = pq.read_table(path, columns=["coord_ra", "coord_dec", "tract"], filters=[("detect_isPrimary", "==", True)])
    return t.to_pandas(ignore_metadata=True)


def _head(path):
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        out = pd.DataFrame({c: np.asarray(d[c]).astype(np.asarray(d[c]).dtype.newbyteorder("=")) for c in HEAD_COLS})
    out["SNID"] = pd.to_numeric(pd.Series(out["SNID"]).astype(str).str.strip()).astype(np.int64)
    out["class"] = re.search(r"ROMAN\+LSST_LARGE_([^/]+)/", path).group(1)
    return out


def read_dat(path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("VARNAMES:"):
                names = line.split()[1:]
                break
    d = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                    names=["tag"] + names, usecols=["tag"] + DAT_COLS, compression="gzip", low_memory=False)
    d = d[d.tag == "SN:"].drop(columns="tag")
    for c in DAT_COLS:
        d[c] = pd.to_numeric(d[c])
    return d.astype({"CID": np.int64, "GALID": np.int64, "GENTYPE": np.int32})


def _pixel(args):
    files, hpix, ids, comp_bands, thr = args
    ids = pa.array(ids, pa.int64())
    g = pq.read_table(files["galaxy"][hpix], columns=GAL_COLS)
    g = g.filter(pc.is_in(g["galaxy_id"], value_set=ids)).to_pandas(ignore_metadata=True)
    f = pq.read_table(files["flux"][hpix], columns=["galaxy_id"] + [f"lsst_flux_{b}" for b in "ugrizy"])
    f = f.filter(pc.is_in(f["galaxy_id"], value_set=ids)).to_pandas(ignore_metadata=True)
    g = g.merge(f, on="galaxy_id", how="left")
    for b in "ugrizy":
        with np.errstate(divide="ignore", invalid="ignore"):
            g[f"host_mag_{b}"] = ZP[b] - 2.5 * np.log10(g.pop(f"lsst_flux_{b}").to_numpy(float))
    frac = {(b, c): np.full(len(g), np.nan) for b in comp_bands for c in ("disk", "bulge", "knot")}
    if len(g):
        import h5py
        with h5py.File(files["sed"][hpix], "r") as h:
            wave = h["meta/wave_list"][:] / 10.0                        # nm
            w = {}
            for b in comp_bands:
                t = np.loadtxt(thr[b])
                w[b] = np.interp(wave, t[:, 0], t[:, 1], left=0, right=0) * np.gradient(wave) / wave
            grp = h["galaxy"]
            for i, gid in enumerate(g.galaxy_id.to_numpy()):
                key = str(gid // 100000)
                if key not in grp or str(gid) not in grp[key]:
                    continue
                sed = grp[key][str(gid)][:].astype(float)               # (disk, bulge, knot) x wave, fnu
                for b in comp_bands:
                    ph = sed @ w[b]
                    tot = ph.sum()
                    if tot > 0:
                        for j, c in enumerate(("disk", "bulge", "knot")):
                            frac[(b, c)][i] = ph[j] / tot
    for (b, c), v in frac.items():
        g[f"host_frac_{c}_{b}"] = v
    return g


def classify(paths, rules):
    out = {k: {} for k in rules}
    for p in paths:
        for k, rx in rules.items():
            m = re.fullmatch(rx, os.path.basename(p))
            if m:
                out[k][int(m.group(1)) if m.groups() else p] = p
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth-index", required=True, help="csv index of the snana parquet, the truth dump and the HEAD files")
    ap.add_argument("--galaxy-index", required=True, help="csv index of the skyCatalog galaxy files and LSST throughputs")
    ap.add_argument("--objects", required=True, help="csv index of the objectTable files")
    ap.add_argument("--coverage-radius", type=float, default=60.0, help="arcsec; nearest primary object closer than this")
    ap.add_argument("--component-bands", default="ri", help="bands for the bulge/disk/knot light fractions")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    t0 = time.time()
    tp = pd.read_csv(a.truth_index)["path"].tolist()
    snana = sorted(p for p in tp if re.fullmatch(r"snana_\d+\.parquet", os.path.basename(p)))
    dat = [p for p in tp if os.path.basename(p) == "OU24_TRANSIENT_TRUTH_HEAD.DAT.gz"]
    heads = sorted(p for p in tp if p.endswith("HEAD.FITS.gz"))
    assert len(dat) == 1 and snana and heads, (len(snana), len(dat), len(heads))

    sn = pd.concat([pq.read_table(f, columns=SNANA_COLS).to_pandas(ignore_metadata=True) for f in snana],
                   ignore_index=True).drop_duplicates("id")
    print(f"{len(sn)} skyCatalog transients from {len(snana)} pixels ({time.time() - t0:.0f}s)", file=sys.stderr)
    with Pool(a.workers) as pool:
        obj = pd.concat(pool.map(_objects, pd.read_csv(a.objects)["path"].tolist()), ignore_index=True)
        tree = cKDTree(unit_vectors(obj.coord_ra.to_numpy(), obj.coord_dec.to_numpy()))
        dist, idx = tree.query(unit_vectors(sn.ra.to_numpy(), sn.dec.to_numpy()), k=1, workers=a.workers)
        sn["nearest_object_arcsec"] = np.rad2deg(2 * np.arcsin(dist / 2)) * 3600
        sn["tract"] = obj.tract.to_numpy()[idx]
        sn = sn[sn.nearest_object_arcsec <= a.coverage_radius].reset_index(drop=True)
        print(f"{len(obj)} primary objects; {len(sn)} transients covered ({time.time() - t0:.0f}s)", file=sys.stderr)
        hd = pd.concat(pool.map(_head, heads), ignore_index=True).drop_duplicates("SNID")
    print(f"{len(hd)} SNANA HEAD rows from {len(heads)} files ({time.time() - t0:.0f}s)", file=sys.stderr)
    dt = read_dat(dat[0])

    df = sn.rename(columns={"id": "sn_id", "ra": "sn_ra", "dec": "sn_dec", "z_CMB": "sn_redshift",
                            "host_id": "skycat_host_id", "host_ra": "skycat_host_ra", "host_dec": "skycat_host_dec"})
    df = df.merge(dt.rename(columns={"CID": "sn_id"}), on="sn_id", how="left")
    hd = hd.rename(columns={"SNID": "sn_id", "HOSTGAL_OBJID": "snana_host_objid", "HOSTGAL_DDLR": "snana_host_ddlr",
                            "HOSTGAL_SNSEP": "snana_host_sep", "HOSTGAL_NMATCH": "snana_host_nmatch",
                            "HOSTGAL_MAG_r": "snana_host_mag_r", "HOSTGAL_LOGMASS": "snana_host_logmass",
                            "SIM_HOSTLIB(LOGMASS_TRUE)": "hostlib_logmass_true"})
    df = df.merge(hd, on="sn_id", how="left")
    df["class"] = df["class"].fillna(df.model_name)
    df["host_id"] = df["GALID"].astype("Int64")
    df["hostless"] = ~(df.host_id.fillna(0) > 0)
    df["skycat_host_missing"] = (df.skycat_host_id <= 0) & ~df.hostless
    for c in ("skycat_host_ra", "skycat_host_dec", "host_sn_sep"):
        df.loc[df.skycat_host_missing, c] = np.nan

    # ---- host truth from the galaxy files
    gp = pd.read_csv(a.galaxy_index)["path"].tolist()
    files = classify(gp, {"galaxy": r"galaxy_(\d+)\.parquet", "flux": r"galaxy_flux_(\d+)\.parquet",
                          "sed": r"galaxy_sed_(\d+)\.hdf5"})
    thr = {re.fullmatch(r"total_(\w)\.dat", os.path.basename(p)).group(1): p for p in gp
           if re.fullmatch(r"total_\w\.dat", os.path.basename(p))}
    ids = np.unique(df.loc[~df.hostless, "host_id"].to_numpy(np.int64))
    hpix = ids // 1_000_000_000
    pix = sorted(int(p) for p in np.unique(hpix))
    missing = [p for p in pix if not all(p in files[k] for k in files)]
    print(f"{len(ids)} host ids in pixels {pix}; pixels without all three files: {missing}", file=sys.stderr)
    jobs = [(files, p, ids[hpix == p], list(a.component_bands), thr) for p in pix if p not in missing]
    with Pool(min(a.workers, len(jobs))) as pool:
        h = pd.concat(pool.map(_pixel, jobs), ignore_index=True)
    h = h.rename(columns=RENAME)
    for c in ("bulge", "disk"):
        gabs = np.hypot(h[f"host_g1_{c}"], h[f"host_g2_{c}"])
        h[f"host_q_{c}"] = (1 - gabs) / (1 + gabs)
        h[f"host_beta_{c}"] = 0.5 * np.degrees(np.arctan2(h[f"host_g2_{c}"], h[f"host_g1_{c}"]))
    h["host_logmass_cat"] = np.log10(h.host_stellar_mass)
    for b in a.component_bands:
        h[f"host_bt_{b}"] = h[f"host_frac_bulge_{b}"]
    h["host_id"] = h.host_id.astype("Int64")
    print(f"{len(h)} of {len(ids)} hosts found in the galaxy files ({time.time() - t0:.0f}s)", file=sys.stderr)
    df = df.merge(h, on="host_id", how="left")
    df["host_missing"] = ~df.hostless & df.host_ra.isna()
    df["sn_id"] = df.sn_id.astype(np.int64)
    for c in ("skycat_host_id", "GALID", "SIM_HOSTLIB_GALID", "snana_host_objid"):
        df[c] = df[c].astype("Int64")
    df = df.sort_values("sn_id").reset_index(drop=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), a.output, compression="zstd")
    print(f"wrote {len(df)} transients, {int(df.hostless.sum())} hostless, {int(df.host_missing.sum())} host not in "
          f"the galaxy files ({time.time() - t0:.0f}s)", file=sys.stderr)
    print(df.groupby("class").agg(n=("sn_id", "size"), z_lt_1=("sn_redshift", lambda z: int((z < 1).sum()))).to_string(),
          file=sys.stderr)


if __name__ == "__main__":
    main()
