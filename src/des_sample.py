"""The DES-SN5YR Dovekie SN Ia sample (IDSURVEY = 10, 1,623 SNe) with its published host matched to DES DR2.

The published host (HOSTGAL_RA / HOSTGAL_DEC, and HOSTGAL2_* for the second-ranked galaxy) is the DES-SN DLR
assignment on the deep SN-field coadds (Wiseman et al. 2020, Qu et al. 2024), visually vetted: a DLR-derived
reference, not truth. It is linked to a DES DR2 object (any object, whatever its star/galaxy class) under the
decision reference_match: pos_1arcsec / pos_3arcsec = the nearest DR2 object within 1" / 3" of the published
position; frac_semimajor = the object, among those within 3", whose r-band half-light ellipse (FLUX_RADIUS /
sqrt(q) along A, B, THETA) puts the published position at the smallest elliptical radius, if that is < 0.5.
host_objectId = -1 when nothing matches, so the published host counts as absent from DR2 (cat1).

Writes sn_id (CID), sn_ra, sn_dec (HEAD file), sn_redshift (zHD), field, the published host's position, DES
magnitudes and d_DLR, host_objectId and host_match_sep_arcsec, the same for HOSTGAL2 (host2_*), host_mag_r
(the published HOSTGAL_MAG_r, for recovery-versus-magnitude plots) and HOST_CONFUSION.
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.table import Table
from scipy.spatial import cKDTree

import des_lib as D

MATCH = {"pos_1arcsec": 1.0, "pos_3arcsec": 3.0, "frac_semimajor": 3.0}
FRAC_SEMIMAJOR = 0.5


def read_head(path):
    cols = ["SNID", "RA", "DEC", "HOSTGAL_NMATCH", "HOSTGAL_OBJID", "HOSTGAL_RA", "HOSTGAL_DEC",
            "HOSTGAL_SNSEP", "HOSTGAL_DDLR", "HOSTGAL_CONFUSION", "HOSTGAL_MAG_g", "HOSTGAL_MAG_r", "HOSTGAL_MAG_i",
            "HOSTGAL_MAG_z", "HOSTGAL_LOGMASS", "HOSTGAL_SPECZ", "HOSTGAL2_OBJID", "HOSTGAL2_RA", "HOSTGAL2_DEC",
            "HOSTGAL2_SNSEP", "HOSTGAL2_DDLR", "HOSTGAL2_MAG_r"]
    h = Table.read(path)[cols].to_pandas()
    for c in ("SNID",):
        h[c] = h[c].apply(lambda s: s.decode() if isinstance(s, bytes) else s).str.strip()
    h["SNID"] = h.SNID.astype(np.int64)
    return h


def match(ref_ra, ref_dec, obj, tree, mode):
    """objectId and separation (arcsec) of the DR2 object matched to each reference position (-1 / NaN if none)."""
    ok = (ref_ra > -90) & (ref_dec > -90)
    oid = np.full(len(ref_ra), -1, np.int64)
    sep = np.full(len(ref_ra), np.nan)
    v = D.unit_vectors(ref_ra[ok], ref_dec[ok])
    r = MATCH[mode]
    hits = tree.query_ball_point(v, D.chord(r))
    ids = obj.objectId.to_numpy()
    ou = D.unit_vectors(obj.coord_ra, obj.coord_dec)
    for i, k in zip(np.flatnonzero(ok), hits):
        if not k:
            continue
        k = np.asarray(k)
        s = D.arcsec_from_chord(np.linalg.norm(ou[k] - D.unit_vectors(ref_ra[i], ref_dec[i]), axis=1))
        if mode == "frac_semimajor":
            from sbw_host.core import sky_to_plane
            c = obj.iloc[k]
            x, y = sky_to_plane(c.coord_ra.to_numpy(), c.coord_dec.to_numpy(), ref_ra[i], ref_dec[i])
            a, b, ellok, phi = D.dlr_ellipse(c, "r", "half_light", None, None, None)
            t = np.deg2rad(phi)
            u, w = -x * np.cos(t) - y * np.sin(t), x * np.sin(t) - y * np.cos(t)
            with np.errstate(invalid="ignore", divide="ignore"):
                rr = np.where(ellok, np.hypot(u / a, w / b), np.inf)
            j = int(np.argmin(rr))
            if rr[j] < FRAC_SEMIMAJOR:
                oid[i], sep[i] = ids[k[j]], s[j]
        else:
            j = int(np.argmin(s))
            oid[i], sep[i] = ids[k[j]], s[j]
    return oid, sep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hosts", required=True, help="directory with DES-Dovekie_Metadata.csv and DES-SN5YR_DES_HEAD.FITS.gz")
    ap.add_argument("--objects", required=True, help="des_dr2_objects cache (data/des_dr2_candidates.parquet)")
    ap.add_argument("--reference-match", required=True, choices=sorted(MATCH))
    ap.add_argument("--record-decisions", action="store_true", help="write the active decisions into the parquet metadata")
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    m = pd.read_csv(f"{a.hosts}/DES-Dovekie_Metadata.csv", sep=r"\s+", comment="#")
    m = m[(m["VARNAMES:"] == "SN:") & (m.IDSURVEY.astype(int) == 10)]
    m = m[["CID", "zHD", "FIELD", "HOST_RA", "HOST_DEC", "HOST_DDLR"]].astype({"CID": np.int64, "zHD": float})
    h = read_head(f"{a.hosts}/DES-SN5YR_DES_HEAD.FITS.gz")
    df = m.merge(h, left_on="CID", right_on="SNID", how="left", validate="1:1")
    if df.SNID.isna().any():
        raise SystemExit("Dovekie SNe missing from the HEAD file")
    obj = D.candidate_columns(pq.read_table(a.objects).to_pandas(ignore_metadata=True))
    tree = cKDTree(D.unit_vectors(obj.coord_ra, obj.coord_dec))
    out = pd.DataFrame({"sn_id": df.CID.to_numpy(np.int64), "sn_ra": df.RA.to_numpy(float), "sn_dec": df.DEC.to_numpy(float),
                        "sn_redshift": df.zHD.to_numpy(float), "field": df.FIELD.to_numpy()})
    for p, pre in (("host", "HOSTGAL"), ("host2", "HOSTGAL2")):
        ra, dec = df[f"{pre}_RA"].to_numpy(float), df[f"{pre}_DEC"].to_numpy(float)
        oid, sep = match(ra, dec, obj, tree, a.reference_match)
        out[f"{p}_ra"], out[f"{p}_dec"] = np.where(ra > -90, ra, np.nan), np.where(dec > -90, dec, np.nan)
        out[f"{p}_objectId"], out[f"{p}_match_sep_arcsec"] = oid, sep
        out[f"{p}_des_objid"] = df[f"{pre}_OBJID"].to_numpy(np.int64)
        out[f"{p}_snsep"] = df[f"{pre}_SNSEP"].to_numpy(float)
        out[f"{p}_ddlr_published"] = df[f"{pre}_DDLR"].to_numpy(float)
        out[f"{p}_mag_r_published"] = df[f"{pre}_MAG_r"].to_numpy(float)
    out["host_mag_r"] = out.host_mag_r_published.where(out.host_mag_r_published.between(10, 35))
    for k in ("HOSTGAL_CONFUSION", "HOSTGAL_NMATCH", "HOSTGAL_LOGMASS", "HOSTGAL_SPECZ"):
        out[k.lower().replace("hostgal_", "host_")] = df[k].to_numpy()
    # how the matched DR2 object is classed (the dominant prototype failure: faint hosts classed as stars)
    ecc = obj.set_index("objectId").extended_class_coadd
    out["host_extended_class_coadd"] = ecc.reindex(out.host_objectId).fillna(-99).to_numpy(np.int64)
    out = out.sort_values("sn_id").reset_index(drop=True)
    out.attrs = {"reference_match": a.reference_match, "n_sn": str(len(out)),
                 "frac_host_matched": str(float((out.host_objectId >= 0).mean())),
                 "frac_host2_present": str(float(out.host2_ra.notna().mean()))}
    if a.record_decisions:
        import json as _json
        import pyarrow as _pa
        import pyarrow.parquet as _pq
        tab = _pa.Table.from_pandas(out, preserve_index=False)
        tab = tab.replace_schema_metadata({**(tab.schema.metadata or {}),
                                           b"decisions": _json.dumps({"reference_match": a.reference_match}).encode()})
        _pq.write_table(tab, a.output)
    else:
        out.to_parquet(a.output, index=False)
    print(f"{len(out)} SNe; published host matched to DR2 ({a.reference_match}) for "
          f"{100 * (out.host_objectId >= 0).mean():.2f}%; median match sep {np.nanmedian(out.host_match_sep_arcsec):.3f}\"; "
          f"HOSTGAL2 present {out.host2_ra.notna().sum()}, matched {(out.host2_objectId >= 0).sum()}", file=sys.stderr)


if __name__ == "__main__":
    main()
