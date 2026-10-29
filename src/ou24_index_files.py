"""Write an index (path, size, mtime) of the OU2024 truth files matching one or more globs, which
astra.yaml declares in place of the files themselves (ou24_sn_truth, ou24_host_truth), like
src/dp02_index_objects.py does for DP0.2.

    python src/ou24_index_files.py [--sha256] <out.csv> <glob> [<glob> ...]

--sha256 adds a content hash column (for small files such as the DP2 DELIGHT tables, data/dp2_delight_yjkang_index.csv).
"""
import csv, glob, hashlib, os, sys
sha = "--sha256" in sys.argv
args = [a for a in sys.argv[1:] if a != "--sha256"]
out, patterns = args[0], args[1:]
paths = sorted({p for pat in patterns for p in glob.glob(pat)})
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["path", "size", "mtime"] + ["sha256"] * sha)
    for p in paths:
        st = os.stat(p)
        w.writerow([p, st.st_size, int(st.st_mtime)] + [hashlib.sha256(open(p, "rb").read()).hexdigest()] * sha)
print(f"{len(paths)} files, {sum(os.path.getsize(p) for p in paths) / 1e9:.2f} GB", file=sys.stderr)
