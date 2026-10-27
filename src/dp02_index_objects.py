"""Write the index of DP0.2 objectTable tract files (path, size, mtime) that
astra.yaml declares as the dp02_object input."""
import glob, os, sys, csv
d = sys.argv[1] if len(sys.argv) > 1 else "/global/cfs/cdirs/lsst/shared/rubin/DP0.2/objectTable"
out = sys.argv[2] if len(sys.argv) > 2 else "data/dp02_objectTable_index.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["path", "size", "mtime"])
    for p in sorted(glob.glob(os.path.join(d, "objectTable_tract_*.parq"))):
        st = os.stat(p); w.writerow([p, st.st_size, int(st.st_mtime)])
