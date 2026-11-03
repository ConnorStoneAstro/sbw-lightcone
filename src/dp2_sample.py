"""The DP2 transient sample chosen by decision dp2_sample, as a table keyed like the DP0.2 sample.

Reads the EDP2 transient working-group file named by the option from the index of those files
(data/dp2_sn_candidates_index.csv), after checking that its size and mtime still match the index, so a
re-release by the collaborators fails loudly instead of silently changing the sample. Writes sn_id (the
DiaObject id), sn_ra, sn_dec, host_objectId = -1 (no truth) and the light-curve fit columns the file has
(SALT3 z, x1, c, mB, FITPROB, ... for rc_salt).

parsnip_ddf_clean and parsnip_ddf_galaxy are parsnip_ddf with stars and AGN removed, from the Legacy Surveys
DR10 south sweeps (--legacy-sweeps; pre-transient imaging) within MATCH_ARCSEC of the transient:

* gaia_star: a source with REF_CAT starting with "G" and Gaia parallax S/N > 3 or proper-motion S/N > 5;
* wise_agn: a TYPE PSF source, not a Gaia star, with W1 and W2 S/N > 3 and W1 - W2 >= 0.8 (Vega; W1_AB - 2.699,
  W2_AB - 3.339);
* legacy_psf_1arcsec: any TYPE PSF source (DUP rows ignored).

parsnip_ddf_clean drops gaia_star or wise_agn; parsnip_ddf_galaxy also drops legacy_psf_1arcsec. Every row
keeps the three flags, gaia_ref_no_astrometry (a Gaia-referenced PSF source that is neither a star by
astrometry nor AGN-coloured; kept, reported), legacy_objectId (BRICKID * 1e6 + OBJID of the source that
flagged the transient, star before AGN before PSF, else the nearest source within the radius; -1 if none),
legacy_type, legacy_sep_arcsec, legacy_n_1arcsec and ddf_field. The parquet attrs key
dp2_sample_accounting (a JSON string) holds the counts removed by each rule, overall and per field.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.io import fits
from scipy.spatial import cKDTree

FILES = {"rc_salt": "RC_sample.parquet", "md_final": "final_sample.parquet", "parsnip_ddf": "parsnip_ddf.parquet",
         "rc_all": "RC_lcfit_input_122k.parquet",
         "parsnip_ddf_clean": "parsnip_ddf.parquet", "parsnip_ddf_galaxy": "parsnip_ddf.parquet"}
CLEAN = {"parsnip_ddf_clean": ("gaia_star", "wise_agn"),
         "parsnip_ddf_galaxy": ("gaia_star", "wise_agn", "legacy_psf_1arcsec")}
MATCH_ARCSEC = 1.0
# DDF centres (deg) for the per-field accounting; each transient goes to the nearest
DDF = {"ELAIS-S1": (9.45, -44.0), "XMM-LSS": (35.57, -4.82), "ECDFS": (53.13, -28.10), "EDFS": (61.0, -48.5),
       "COSMOS": (150.11, 2.18)}
LEG = ["BRICKID", "OBJID", "TYPE", "RA", "DEC", "REF_CAT", "PARALLAX", "PARALLAX_IVAR", "PMRA", "PMRA_IVAR", "PMDEC",
       "PMDEC_IVAR", "FLUX_W1", "FLUX_W2", "FLUX_IVAR_W1", "FLUX_IVAR_W2"]
# (id, ra, dec) column names per file
POS = {"rc_salt": ("diaObjectId", "ra", "dec"),
       "md_final": ("diaObjectId_dia_object_lc", "ra_dia_object_lc", "dec_dia_object_lc"),
       "parsnip_ddf": ("diaObjectId", "ra", "dec"),
       "rc_all": ("diaObjectId", "ra", "dec")}
KEEP = ["z", "z_err", "PKMJD", "x1", "x1ERR", "c", "cERR", "mB", "mBERR", "FITPROB", "MWEBV", "SNRMAX1", "SNRSUM",
        "NOBS_DETECT", "N_NIGHTS_SNR5", "N_BANDS_SNR5", "in_MD", "MD_set", "MD_parsnip_p_SNIa",
        "predicted_redshift", "count"]


def unit_vectors(ra, dec):
    ra, dec = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def ddf_field(ra, dec):
    names = list(DDF)
    c = unit_vectors([DDF[n][0] for n in names], [DDF[n][1] for n in names])
    return np.array(names)[np.argmax(unit_vectors(ra, dec) @ c.T, axis=1)]


def legacy_within(sn, sweeps_csv, radius):
    """Legacy DR10 sources (not DUP) within radius arcsec of each transient, one row per pair, with the
    star / AGN / PSF classification."""
    from legacy_index_sweeps import tiles_near
    paths = pd.read_csv(sweeps_csv)["path"].tolist()
    need = set().union(*[tiles_near(r, d, radius / 3600.0) for r, d in zip(sn.sn_ra, sn.sn_dec)])
    missing = need - {os.path.basename(p) for p in paths}
    if missing:
        raise SystemExit(f"sweep tiles {sorted(missing)} are needed but not in {sweeps_csv}")
    xyz = unit_vectors(sn.sn_ra, sn.sn_dec)
    chord = 2 * np.sin(np.deg2rad(radius / 3600.0) / 2)
    parts = []
    for p in paths:
        if os.path.basename(p) not in need:
            continue
        with fits.open(p, memmap=True) as h:
            t = h[1].data
            hits = cKDTree(unit_vectors(t["RA"], t["DEC"])).query_ball_point(xyz, chord)
            pairs = [(i, j) for i, js in enumerate(hits) for j in js]
            if not pairs:
                continue
            si, oj = np.array(pairs).T
            native = lambda x: np.asarray(x).astype(np.asarray(x).dtype.newbyteorder("="))   # FITS is big-endian
            c = pd.DataFrame({k: native(t[k][oj]) for k in LEG})
        c.insert(0, "sn_id", sn.sn_id.to_numpy()[si])
        c["sep_arcsec"] = np.rad2deg(2 * np.arcsin(np.linalg.norm(xyz[si] - unit_vectors(c.RA, c.DEC), axis=1) / 2)) * 3600
        parts.append(c)
    c = pd.concat(parts, ignore_index=True)
    c["TYPE"] = c.TYPE.astype(str).str.strip()
    c = c[c.TYPE != "DUP"].drop_duplicates(["sn_id", "BRICKID", "OBJID"]).reset_index(drop=True)
    c["legacy_objectId"] = c.BRICKID.astype(np.int64) * 1_000_000 + c.OBJID.astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        plx = c.PARALLAX * np.sqrt(c.PARALLAX_IVAR)
        pm = np.hypot(c.PMRA * np.sqrt(c.PMRA_IVAR), c.PMDEC * np.sqrt(c.PMDEC_IVAR))
        w12 = -2.5 * np.log10(c.FLUX_W1 / c.FLUX_W2) + (3.339 - 2.699)
        w_ok = (c.FLUX_W1 * np.sqrt(c.FLUX_IVAR_W1) > 3) & (c.FLUX_W2 * np.sqrt(c.FLUX_IVAR_W2) > 3)
    gaia = c.REF_CAT.astype(str).str.strip().str.startswith("G")
    psf = c.TYPE == "PSF"
    c["gaia_star"] = gaia & ((plx > 3) | (pm > 5))
    c["wise_agn"] = psf & ~c.gaia_star & w_ok & (w12 >= 0.8)
    c["legacy_psf_1arcsec"] = psf
    c["gaia_ref_no_astrometry"] = gaia & psf & ~c.gaia_star & ~c.wise_agn
    # which source names the transient: star, then AGN, then PSF, then the nearest
    c["rank"] = np.select([c.gaia_star, c.wise_agn, psf], [0, 1, 2], 3)
    return c.sort_values(["sn_id", "rank", "sep_arcsec"], kind="stable")


def classify(df, sweeps_csv, sample):
    c = legacy_within(df, sweeps_csv, MATCH_ARCSEC)
    flags = ["gaia_star", "wise_agn", "legacy_psf_1arcsec", "gaia_ref_no_astrometry"]
    g = c.groupby("sn_id")
    any_ = g[flags].any()
    first = g.head(1).set_index("sn_id")
    for f in flags:
        df[f] = any_[f].reindex(df.sn_id).fillna(False).astype(bool).to_numpy()
    df["legacy_objectId"] = first.legacy_objectId.reindex(df.sn_id).fillna(-1).astype(np.int64).to_numpy()
    df["legacy_type"] = first.TYPE.reindex(df.sn_id).fillna("none").to_numpy()
    df["legacy_sep_arcsec"] = first.sep_arcsec.reindex(df.sn_id).to_numpy()
    df["legacy_n_1arcsec"] = g.size().reindex(df.sn_id).fillna(0).astype(int).to_numpy()
    df["ddf_field"] = ddf_field(df.sn_ra, df.sn_dec)
    # sequential accounting: each rule counts the transients it removes after the previous ones
    acc = {}
    for name, sub in [("all", df)] + [(f, df[df.ddf_field == f]) for f in DDF if (df.ddf_field == f).any()]:
        keep = np.ones(len(sub), bool)
        row = {"parsnip_ddf": int(len(sub))}
        for rule in ("gaia_star", "wise_agn", "legacy_psf_1arcsec"):
            rm = keep & sub[rule].to_numpy()
            row[f"removed_{rule}"] = int(rm.sum())
            keep &= ~rm
            row[f"after_{rule}"] = int(keep.sum())
        clean = ~(sub.gaia_star | sub.wise_agn).to_numpy()
        row["gaia_ref_no_astrometry_kept_in_clean"] = int((clean & sub.gaia_ref_no_astrometry.to_numpy()).sum())
        row["no_legacy_source_1arcsec"] = int((sub.legacy_n_1arcsec == 0).sum())
        acc[name] = row
    drop = np.zeros(len(df), bool)
    for f in CLEAN[sample]:
        drop |= df[f].to_numpy()
    return df[~drop].reset_index(drop=True), acc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", required=True, help="csv (path, size, mtime) of the EDP2 sample files")
    ap.add_argument("--sample", required=True, choices=sorted(FILES))
    ap.add_argument("--legacy-sweeps", help="legacy_dr10_sweeps index csv; read only for the parsnip_ddf_* "
                    "star/AGN-cleaned options")
    ap.add_argument("--record-decisions", action="store_true", help="write the active decisions into the parquet metadata")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    idx = pd.read_csv(a.index)
    row = idx[idx.path.str.endswith("/" + FILES[a.sample])]
    if len(row) != 1:
        raise SystemExit(f"{FILES[a.sample]} is not in {a.index}")
    path, size, mtime = row.iloc[0][["path", "size", "mtime"]]
    st = os.stat(path)
    if st.st_size != size or int(st.st_mtime) != mtime:
        raise SystemExit(f"{path} changed since the index was written (size {st.st_size} vs {size}, "
                         f"mtime {int(st.st_mtime)} vs {mtime}); rewrite the index deliberately")
    sid, ra, dec = POS[a.sample.replace("_clean", "").replace("_galaxy", "")]
    names = pq.read_schema(path).names
    df = pq.read_table(path, columns=[sid, ra, dec] + [c for c in KEEP if c in names]).to_pandas(ignore_metadata=True)
    df = df.rename(columns={sid: "sn_id", ra: "sn_ra", dec: "sn_dec"})
    df["sn_id"] = df["sn_id"].astype(np.int64)
    if df.sn_id.duplicated().any():
        raise SystemExit(f"duplicate DiaObject ids in {path}")
    df["host_objectId"] = np.int64(-1)
    df = df.sort_values("sn_id").reset_index(drop=True)
    meta = {"dp2_sample": a.sample, "source": path}
    if a.record_decisions:
        meta["decisions"] = json.dumps({"dp2_sample": a.sample})
    if a.sample in CLEAN:
        if not a.legacy_sweeps:
            raise SystemExit(f"{a.sample} needs --legacy-sweeps")
        df, acc = classify(df, a.legacy_sweeps, a.sample)
        meta["dp2_sample_accounting"] = json.dumps(acc)
        print(json.dumps(acc, indent=1), file=sys.stderr)
    df.attrs = meta
    tab = pa.Table.from_pandas(df, preserve_index=False)
    tab = tab.replace_schema_metadata({**(tab.schema.metadata or {}), **{k.encode(): str(v).encode() for k, v in meta.items()}})
    pq.write_table(tab, a.output)
    print(f"{a.sample}: {len(df)} transients from {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
