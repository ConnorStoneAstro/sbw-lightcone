"""Ten random positions per DES SN, uniform in area in the 45-90 arcsec annulus around it, so they share the
SN's deep-field depth and their 30" search cone stays inside the SN's 2' des_dr2_objects cone. Ids are
-(100 * k + j + 1) for the k-th SN (by sn_id) and its j-th position; parent_sn_id records the SN. Seeded."""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

SEED = 20261029
R_IN, R_OUT = 45.0, 90.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--per-sn", type=int, default=10)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    sn = pd.read_parquet(a.sample, columns=["sn_id", "sn_ra", "sn_dec"]).sort_values("sn_id").reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    n, m = len(sn), a.per_sn
    r = np.sqrt(rng.uniform(R_IN ** 2, R_OUT ** 2, (n, m))) / 3600.0
    th = rng.uniform(0, 2 * np.pi, (n, m))
    ra0, dec0 = np.deg2rad(sn.sn_ra.to_numpy())[:, None], np.deg2rad(sn.sn_dec.to_numpy())[:, None]
    rr = np.deg2rad(r)
    # destination point at angular distance rr and bearing th (from North through East)
    dec = np.arcsin(np.sin(dec0) * np.cos(rr) + np.cos(dec0) * np.sin(rr) * np.cos(th))
    ra = ra0 + np.arctan2(np.sin(th) * np.sin(rr) * np.cos(dec0), np.cos(rr) - np.sin(dec0) * np.sin(dec))
    k, j = np.meshgrid(np.arange(n), np.arange(m), indexing="ij")
    out = pd.DataFrame({"sn_id": -(100 * k + j + 1).ravel().astype(np.int64),
                        "sn_ra": np.rad2deg(ra).ravel() % 360.0, "sn_dec": np.rad2deg(dec).ravel(),
                        "parent_sn_id": np.repeat(sn.sn_id.to_numpy(np.int64), m),
                        "offset_arcsec": (r * 3600.0).ravel()})
    out["host_objectId"] = np.int64(-1)
    out.to_parquet(a.output, index=False)


if __name__ == "__main__":
    main()
