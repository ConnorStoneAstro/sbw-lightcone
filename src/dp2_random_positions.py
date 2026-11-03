"""Random sky positions for the DP2 hostless control: ten per sample transient, uniform in area in the
45-90 arcsec annulus around it (dp02_lib.annulus_positions, the design every tier shares since 2 Nov 2026),
so they share the transient's field (DDF or WFD) and depth. Ids are -(100 k + j + 1) for the k-th transient
(by sn_id) and its j-th position; parent_sn_id records the transient. Before, the positions were drawn
uniformly over the transient's own patch."""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from dp02_lib import annulus_positions


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True, help="dp2_sample")
    ap.add_argument("--per-sn", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    if a.per_sn > 99:
        raise SystemExit("--per-sn must be < 100 (id scheme)")
    sn = pd.read_parquet(a.sample, columns=["sn_id", "sn_ra", "sn_dec"])
    out = annulus_positions(sn, a.per_sn, a.seed)
    out.to_parquet(a.output, index=False)
    print(f"wrote {len(out)} random positions for {len(sn)} transients", file=sys.stderr)


if __name__ == "__main__":
    main()
