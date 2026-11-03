"""Random sky positions for the simulated hostless control: ten per sample SN, uniform in area in the
45-90 arcsec annulus around it (dp02_lib.annulus_positions, the design every tier shares since 2 Nov 2026).
They share the SN's field and depth; the 30 arcsec search cone around a position never contains the SN.
Ids are -(100 k + j + 1) for the k-th SN (by sn_id) and its j-th position; parent_sn_id records the SN.
Before, one position per truth SN was drawn uniformly over the SN's tract."""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from dp02_lib import annulus_positions


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, help="dp02_sample or ou24_sample (sn_id, sn_ra, sn_dec)")
    ap.add_argument("--per-sn", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    if a.per_sn > 99:
        raise SystemExit("--per-sn must be < 100 (id scheme)")
    sn = pd.read_parquet(a.sample, columns=["sn_id", "sn_ra", "sn_dec"])
    out = annulus_positions(sn, a.per_sn, a.seed)
    out.to_parquet(a.output, index=False, compression="zstd")
    print(f"wrote {len(out)} random positions for {len(sn)} SNe", file=sys.stderr)


if __name__ == "__main__":
    main()
