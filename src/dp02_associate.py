"""SBW and DLR association of every sample SN against its DP0.2 candidates."""
from __future__ import annotations

import argparse

import pyarrow as pa
import pyarrow.parquet as pq

from dp02_lib import DLR_MAX, PROFILES, RADIUS, associate_all, load_candidates


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--profile", required=True, choices=PROFILES + ("stamps",))
    ap.add_argument("--band", required=True, choices=list("ugrizy"))
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--dlr-max", required=True, choices=sorted(DLR_MAX))
    ap.add_argument("--dlr-radius", required=True, choices=["half_light", "moments"])
    ap.add_argument("--extendedness", required=True, choices=["ext05", "none"])
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    if a.profile == "stamps":
        raise SystemExit("sbw_profile=stamps needs deblended stamps, which do not exist for DP0.2 on NERSC")
    sample = pq.read_table(a.sample).to_pandas()
    cand = load_candidates(a.candidates, bands=a.band)
    out = associate_all(sample, cand, band=a.band, profile=a.profile, search_radius=a.search_radius,
                        dlr_max=a.dlr_max, dlr_radius=a.dlr_radius, extendedness=a.extendedness, workers=a.workers)
    tab = pa.Table.from_pandas(out, preserve_index=False)
    tab = tab.replace_schema_metadata({k.encode(): str(v).encode() for k, v in out.attrs.items()})
    pq.write_table(tab, a.output, compression="zstd")


if __name__ == "__main__":
    main()
