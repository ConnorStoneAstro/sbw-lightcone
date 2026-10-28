"""Write an index (path, size, mtime) of the OU2024 truth files matching one or more globs, which
astra.yaml declares in place of the files themselves (ou24_sn_truth, ou24_host_truth), like
src/dp02_index_objects.py does for DP0.2.

    python src/ou24_index_files.py <out.csv> <glob> [<glob> ...]
"""
import csv, glob, os, sys
out, patterns = sys.argv[1], sys.argv[2:]
paths = sorted({p for pat in patterns for p in glob.glob(pat)})
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["path", "size", "mtime"])
    for p in paths:
        st = os.stat(p); w.writerow([p, st.st_size, int(st.st_mtime)])
print(f"{len(paths)} files, {sum(os.path.getsize(p) for p in paths) / 1e9:.2f} GB", file=sys.stderr)
