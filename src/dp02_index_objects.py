"""Write an index (path, size, mtime) of a large read-only file collection, which
astra.yaml declares in place of the files themselves: the DP0.2 objectTable tract
files (dp02_object) and the cosmoDC2 image healpix files (cosmodc2_hosts).

    python src/dp02_index_objects.py [dir] [out] [glob]
"""
import glob, os, sys, csv
d = sys.argv[1] if len(sys.argv) > 1 else "/global/cfs/cdirs/lsst/shared/rubin/DP0.2/objectTable"
out = sys.argv[2] if len(sys.argv) > 2 else "data/dp02_objectTable_index.csv"
pattern = sys.argv[3] if len(sys.argv) > 3 else "objectTable_tract_*.parq"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["path", "size", "mtime"])
    for p in sorted(glob.glob(os.path.join(d, pattern))):
        st = os.stat(p); w.writerow([p, st.st_size, int(st.st_mtime)])
