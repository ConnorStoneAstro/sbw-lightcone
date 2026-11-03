"""SBW and DLR association of every sample SN against its DP0.2 candidates (OU2024 the same; DP2 with
``--tier dp2``, which swaps in the MultiProFit light models of src/dp2_lib.py)."""
from __future__ import annotations

import argparse
import json

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

import dp02_lib as L
from dp02_lib import DLR_MAX, PROFILES, QUALITY, RADIUS, associate_all


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sample", help="dp02_sample table (SNe with their true host objectId)")
    src.add_argument("--positions", help="positions with no true host (sn_id, sn_ra, sn_dec), e.g. dp02_random_positions")
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--profile", required=True, choices=PROFILES + ("stamps",))
    ap.add_argument("--band", required=True, choices=list("ugrizy"))
    ap.add_argument("--search-radius", required=True, choices=sorted(RADIUS))
    ap.add_argument("--dlr-max", required=True, choices=sorted(DLR_MAX))
    ap.add_argument("--dlr-radius", required=True, choices=["half_light", "moments"])
    ap.add_argument("--extendedness", required=True, choices=["ext05", "none"])
    ap.add_argument("--candidate-quality", required=True, choices=QUALITY)
    ap.add_argument("--tier", choices=["dp02", "ou24", "dp2"], default="dp02")
    ap.add_argument("--exclude-suspect", choices=["keep", "drop"], default="keep", help="DP2 only (sbw_host suspect fits)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--output", required=True)
    a = ap.parse_args(argv)
    if a.profile == "stamps":
        raise SystemExit("sbw_profile=stamps needs deblended stamps, which do not exist on NERSC")
    if a.tier == "dp2":
        import dp2_lib
        dp2_lib.install(a.exclude_suspect)
    recorded = {}
    for path in (a.sample, a.candidates):   # decisions the upstream recipes recorded (--record-decisions)
        md = pq.read_schema(path).metadata or {} if path else {}
        if b"decisions" in md:
            recorded.update(json.loads(md[b"decisions"]))
    if a.sample:
        sample = pq.read_table(a.sample).to_pandas(ignore_metadata=True)
    else:
        sample = pq.read_table(a.positions).to_pandas(ignore_metadata=True)
        sample["host_objectId"] = np.int64(-1)
    cand = L.load_candidates(a.candidates, bands=a.band)
    out = associate_all(sample, cand, band=a.band, profile=a.profile, search_radius=a.search_radius,
                        dlr_max=a.dlr_max, dlr_radius=a.dlr_radius, extendedness=a.extendedness,
                        quality=a.candidate_quality, workers=a.workers)
    if a.tier == "dp2":
        out.attrs["exclude_suspect"] = a.exclude_suspect
    out.attrs.update(recorded)   # so every settings block names the sample it stands on
    tab = pa.Table.from_pandas(out, preserve_index=False)
    tab = tab.replace_schema_metadata({k.encode(): str(v).encode() for k, v in out.attrs.items()})
    pq.write_table(tab, a.output, compression="zstd")


if __name__ == "__main__":
    main()
