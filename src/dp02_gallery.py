"""Gallery of Rubin DP0.2 (or, with ``--tier dp2``, DP2) deep-coadd cutouts with host-candidate markers.

Reads a csv of cases and renders one panel per case: a ``deepCoadd_calexp``
cutout centred on the supernova, north up and east left, with a marker for
every ``<name>_ra`` / ``<name>_dec`` column pair (SBW pick, DLR pick, truth
matches, ...).

Image access (off-RSP, e.g. from a laptop)
------------------------------------------
1. ``dp02_dc2_catalogs.ObsCore`` over TAP (``read:tap`` scope) gives one row per
   patch and band (``dataproduct_subtype = 'lsst.deepCoadd_calexp'``); the
   butler dataset ID is the ``ID`` query parameter of ``access_url``.  The patch
   table for a band is fetched once and cached.
2. The SODA service ``/api/cutout/sync`` with ``ID`` and a square ``POLYGON``
   returns a FITS cutout (``read:image`` scope; a ``read:tap``-only token gets 403).

Tokens are taken from ``--token-file`` / ``--tap-token-file`` /
``--image-token-file``, else ``RSP_TOKEN`` / ``RSP_IMAGE_TOKEN``, else the files
``.rsp_image_tap_token`` (both scopes), ``.rsp_token`` and ``.rsp_image_token``
found in the working directory or any parent.

Usage::

    python gallery.py cases.csv --output gallery.pdf --band r --size 30 --ncol 4

Input columns: ``case_id``, ``sn_ra``, ``sn_dec``, optional ``label``, and any
number of ``<name>_ra`` / ``<name>_dec`` marker pairs (blank = not drawn).
A marker with ``<name>_a``, ``<name>_b`` and ``<name>_phi`` columns (semi-axes in arcsec, angle in degrees
counter-clockwise from +x = West toward North) also gets that ellipse drawn about its position.
With ``--tier dp2`` the patch table is ``ivoa.ObsCore`` rows of ``obs_collection = 'LSST.DP2'`` and subtype
``lsst.deep_coadd``.
With ``--lookup sbw,dlr`` a marker may instead be given as ``<name>_objectId``
(Int64; blank or negative = not drawn) and its position is read from
``dp02_dc2_catalogs.Object`` over TAP.
"""

from __future__ import annotations

import argparse
import io
import os
import random
import re
import textwrap
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from astropy.io import fits  # noqa: E402
from astropy.wcs import WCS  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

RSP_BASE = "https://data.lsst.cloud"
TAP_URL = f"{RSP_BASE}/api/tap"
CUTOUT_URL = f"{RSP_BASE}/api/cutout/sync"
# tier -> (ObsCore table, extra WHERE clause, coadd subtype, title)
TIERS = {"dp02": ("dp02_dc2_catalogs.ObsCore", "", "lsst.deepCoadd_calexp", "DP0.2 deepCoadd"),
         "dp2": ("ivoa.ObsCore", " AND obs_collection = 'LSST.DP2'", "lsst.deep_coadd", "DP2 deep_coadd")}
OBSCORE_TABLE, OBS_WHERE, COADD_SUBTYPE, TITLE = TIERS["dp02"]

SBW_COLOR, DLR_COLOR, GREY = "#2a6fdb", "#e07b39", "#555555"
SN_COLOR = "#d62728"
# name -> (marker, colour, size, legend label).  Open markers of different
# sizes, so coincident picks stay visible as concentric symbols.
MARKER_STYLES = {
    "sbw": ("o", SBW_COLOR, 13, "SBW most probable"),
    "dlr": ("s", DLR_COLOR, 10, "DLR"),
    "truth_match": ("D", "black", 16, "MatchesTruth"),
    "truth_nearest": ("h", GREY, 21, "Nearest to truth host"),
    "truth_pos": ("+", "#999999", 7, "Truth host position"),
    "flip": ("^", "#2ca02c", 11, "SBW most probable, other band"),
    "delight": ("*", "#9467bd", 16, "DELIGHT (CNN position)"),
    "sex": ("v", "#2ca02c", 11, "DELIGHT-SEx host"),
}
EXTRA_STYLES = [("v", "#2ca02c"), ("h", "#9467bd"), ("<", "#8c564b"), (">", "#e377c2"), ("*", "#17becf")]


# --------------------------------------------------------------------- tokens

def _find_up(name: str) -> Path | None:
    for d in [Path.cwd(), *Path.cwd().parents]:
        p = d / name
        if p.is_file():
            return p
    return None


def _read_token(path) -> str:
    text = Path(path).read_text().strip()
    if not text:
        raise RuntimeError(f"Token file {path} is empty.")
    return text


def load_tokens(token_file=None, tap_token_file=None, image_token_file=None) -> tuple[str, str]:
    """Return ``(tap_token, image_token)``; see the module docstring for the search order."""
    both = _read_token(token_file) if token_file else None
    tap = (_read_token(tap_token_file) if tap_token_file else None) or both or os.environ.get("RSP_TOKEN")
    img = (_read_token(image_token_file) if image_token_file else None) or both or os.environ.get("RSP_IMAGE_TOKEN")
    combined = _find_up(".rsp_image_tap_token")
    if not tap:
        p = combined or _find_up(".rsp_token")
        tap = _read_token(p) if p else None
    if not img:
        p = combined or _find_up(".rsp_image_token")
        img = _read_token(p) if p else None
    if not tap or not img:
        raise RuntimeError(
            "Need RSP tokens with read:tap (ObsCore lookup) and read:image (cutouts). Create one at "
            f"{RSP_BASE} (Account -> Security tokens) and save it as .rsp_image_tap_token, or pass "
            "--token-file / --tap-token-file / --image-token-file."
        )
    return tap.strip(), img.strip()


def _session(token: str) -> requests.Session:
    s = requests.Session()
    s.headers["Authorization"] = f"Bearer {token}"
    adapter = requests.adapters.HTTPAdapter(pool_connections=8, pool_maxsize=16)
    s.mount("https://", adapter)
    return s


# -------------------------------------------------------------- patch lookup

_POLY_RE = re.compile(r"POLYGON\s+(?:ICRS\s+)?(.*)", re.I)


def patch_table(band: str, tap_token: str, cache_dir: Path) -> pd.DataFrame:
    """All DP0.2 deep-coadd patches in ``band`` with their ra/dec bounding boxes (cached)."""
    cache = cache_dir / f"obscore_{OBSCORE_TABLE.split('.')[0]}_{COADD_SUBTYPE.split('.')[-1]}_{band}.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    import pyvo

    tap = pyvo.dal.TAPService(TAP_URL, session=_session(tap_token))
    q = (f"SELECT access_url, lsst_tract, lsst_patch, s_ra, s_dec, s_region FROM {OBSCORE_TABLE} "
         f"WHERE dataproduct_subtype = '{COADD_SUBTYPE}' AND lsst_band = '{band}'{OBS_WHERE}")
    df = tap.search(q, maxrec=200000).to_table().to_pandas()
    if df.empty:
        raise RuntimeError(f"ObsCore returned no {COADD_SUBTYPE} rows for band {band!r}.")
    df["dataset_id"] = [parse_qs(urlparse(str(u)).query)["ID"][0] for u in df["access_url"]]
    corners = [np.array(_POLY_RE.match(str(r)).group(1).split(), float).reshape(-1, 2) for r in df["s_region"]]
    for c in corners:   # a patch straddling RA = 0: unwrap to a contiguous range (ra_min < 0)
        if np.ptp(c[:, 0]) > 180:
            c[:, 0] = np.where(c[:, 0] > 180, c[:, 0] - 360, c[:, 0])
    df["ra_min"] = [c[:, 0].min() for c in corners]
    df["ra_max"] = [c[:, 0].max() for c in corners]
    df["dec_min"] = [c[:, 1].min() for c in corners]
    df["dec_max"] = [c[:, 1].max() for c in corners]
    df = df.drop(columns=["s_region"])
    cache_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache)
    return df


def best_patch(patches: pd.DataFrame, ra: float, dec: float) -> pd.Series:
    """Patch containing (ra, dec) with the largest margin to its edges (least clipping)."""
    cosd = np.cos(np.radians(dec))
    margin = np.maximum.reduce([np.minimum.reduce([
        (r - patches["ra_min"].to_numpy()) * cosd, (patches["ra_max"].to_numpy() - r) * cosd,
        dec - patches["dec_min"].to_numpy(), patches["dec_max"].to_numpy() - dec,
    ]) for r in (ra, ra - 360)])
    i = int(np.argmax(margin))
    if margin[i] <= 0:
        raise LookupError(f"No {TITLE} covers ({ra:.5f}, {dec:.5f}).")
    return patches.iloc[i]


def lookup_positions(cases: pd.DataFrame, names: list[str], tap_token: str) -> pd.DataFrame:
    """Fill ``<name>_ra/_dec`` from ``<name>_objectId`` via the Object table (chunked IN queries)."""
    import pyvo

    ids = {}
    for n in names:
        col = cases[f"{n}_objectId"].astype("Int64")
        cases[f"{n}_objectId"] = col
        ids.update({int(v): None for v in col.dropna() if v > 0})
    if not ids:
        return cases
    tap = pyvo.dal.TAPService(TAP_URL, session=_session(tap_token))
    keys = list(ids)
    parts = [tap.search(f"SELECT objectId, coord_ra, coord_dec FROM dp02_dc2_catalogs.Object "
                        f"WHERE objectId IN ({','.join(map(str, keys[i:i + 500]))})").to_table().to_pandas()
             for i in range(0, len(keys), 500)]
    pos = pd.concat(parts).set_index("objectId")
    for n in names:
        oid = cases[f"{n}_objectId"].where(cases[f"{n}_objectId"] > 0)
        cases[f"{n}_ra"] = oid.map(pos["coord_ra"]).astype(float)
        cases[f"{n}_dec"] = oid.map(pos["coord_dec"]).astype(float)
    return cases


# ------------------------------------------------------------------- cutouts

class RateLimiter:
    """At most one request start every ``interval`` seconds, across threads."""

    def __init__(self, interval: float):
        self.interval, self._next, self._lock = interval, 0.0, threading.Lock()

    def wait(self):
        with self._lock:
            now = time.monotonic()
            t = max(now, self._next)
            self._next = t + self.interval
        time.sleep(max(0.0, t - now))


def cache_path(cache_dir: Path, ra: float, dec: float, band: str, size: float) -> Path:
    return cache_dir / f"dp02_{band}_{ra:.6f}_{dec:+.6f}_{size:g}as.fits"


def fetch_cutout(ra, dec, band, size, patches, session, cache_dir, limiter, retries=4):
    """FITS bytes of a ``size``-arcsec square cutout centred on (ra, dec), cached on disk."""
    path = cache_path(cache_dir, ra, dec, band, size)
    if path.exists():
        return path.read_bytes()
    patch = best_patch(patches, ra, dec)
    h = size / 2.0 / 3600.0
    cosd = np.cos(np.radians(dec))
    poly = " ".join(f"{ra + sx * h / cosd:.8f} {dec + sy * h:.8f}"
                    for sx, sy in [(-1, -1), (1, -1), (1, 1), (-1, 1)])
    params = {"ID": patch["dataset_id"], "POLYGON": poly}
    err = None
    for attempt in range(retries):
        limiter.wait()
        try:
            r = session.get(CUTOUT_URL, params=params, timeout=180)
            if r.status_code in (401, 403):
                raise PermissionError(f"cutout service HTTP {r.status_code}: token lacks read:image scope?")
            if r.status_code == 429 or r.status_code >= 500:
                raise requests.HTTPError(f"HTTP {r.status_code}")
            r.raise_for_status()
            if not r.content.startswith(b"SIMPLE"):
                raise ValueError(f"not FITS: {r.content[:80]!r}")
            cache_dir.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(f".tmp{threading.get_ident()}")
            tmp.write_bytes(r.content)
            tmp.replace(path)
            return r.content
        except PermissionError:
            raise
        except (requests.RequestException, ValueError) as e:
            err = e
            time.sleep(2.0 * 2**attempt + random.random())
    raise RuntimeError(f"cutout failed after {retries} tries: {err}")


def read_image(content: bytes) -> tuple[np.ndarray, WCS]:
    with fits.open(io.BytesIO(content)) as hdul:
        for hdu in hdul:
            if hdu.data is not None and np.ndim(hdu.data) == 2:
                return np.asarray(hdu.data, float), WCS(hdu.header).celestial
    raise ValueError("cutout has no 2-D image HDU")


# ------------------------------------------------------------------ plotting

def tangent_offsets(ra, dec, ra0, dec0):
    """Gnomonic offsets (east, north) in arcsec of (ra, dec) about (ra0, dec0)."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi) * 3600.0, np.degrees(eta) * 3600.0


def stretch(data: np.ndarray, kind: str, pmin: float, pmax: float) -> np.ndarray:
    """Clip to the [pmin, pmax] percentiles and map to [0, 1].

    For asinh the softening (linear-to-log knee) sits at 3 sigma of the
    background (MAD estimate), so sky noise stays near the background shade in empty fields while
    bright galaxies are still compressed.
    """
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return np.zeros_like(data)
    lo, hi = np.percentile(finite, [pmin, pmax])
    span = max(hi - lo, 1e-30)
    x = np.clip((data - lo) / span, 0, 1)
    if kind == "asinh":
        sigma = 1.4826 * np.median(np.abs(finite - np.median(finite)))
        soften = float(np.clip(3 * sigma / span, 1e-3, 1.0))
        return np.arcsinh(x / soften) / np.arcsinh(1 / soften)
    if kind == "sqrt":
        return np.sqrt(x)
    return x


def marker_names(columns) -> list[str]:
    cols = set(columns)
    return [c[:-3] for c in columns if c.endswith("_ra") and c != "sn_ra" and c[:-3] + "_dec" in cols]


def marker_style(name: str, i_extra: int):
    if name in MARKER_STYLES:
        return MARKER_STYLES[name]
    m, c = EXTRA_STYLES[i_extra % len(EXTRA_STYLES)]
    return m, c, 13, name


def draw_panel(ax, case, image, error, names, styles, size, args):
    half = size / 2.0
    halo = [patheffects.withStroke(linewidth=2, foreground="black" if args.dark else "white")]
    sn_ra, sn_dec = float(case["sn_ra"]), float(case["sn_dec"])
    if image is not None:
        data, wcs = image
        ny, nx = data.shape
        # pixel-corner grid -> (east, north) offsets about the SN: exact for any WCS orientation
        yy, xx = np.mgrid[-0.5:ny, -0.5:nx]
        cra, cdec = wcs.wcs_pix2world(xx, yy, 0)
        ex, no = tangent_offsets(cra, cdec, sn_ra, sn_dec)
        ax.pcolormesh(ex, no, stretch(data, args.stretch, args.pmin, args.pmax),
                      cmap="gray" if args.dark else "gray_r", shading="flat", rasterized=True)
    else:
        ax.set_facecolor("#eeeeee")
        ax.text(0, 0, textwrap.fill(f"cutout failed: {error}", 34), ha="center", va="center",
                fontsize=6, color="#aa0000")
        ax.set_xlim(half, -half), ax.set_ylim(-half, half), ax.set_aspect("equal")
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_title(_title(case), fontsize=7)
        return
    for n in names:
        r, d = case.get(f"{n}_ra"), case.get(f"{n}_dec")
        if pd.isna(r) or pd.isna(d):
            continue
        x, y = tangent_offsets(float(r), float(d), sn_ra, sn_dec)
        m, c, s, _ = styles[n]
        ea, eb, ephi = (case.get(f"{n}_{k}") for k in ("a", "b", "phi"))
        if not any(pd.isna(v) for v in (ea, eb, ephi)):   # angle from West toward North -> from East (+x here)
            k = 2 * args.ellipse_scale
            ax.add_patch(Ellipse((x, y), k * float(ea), k * float(eb), angle=180.0 - float(ephi), fill=False,
                                 ec=c, lw=1.0, ls="--" if n == "dlr" else "-", zorder=3))
        far = max(abs(x), abs(y))
        if far > 0.95 * half:  # off-panel: arrow at the edge pointing to it, with its distance
            f = 0.93 * half / far
            ax.annotate("", (x * f, y * f), (x * f * 0.8, y * f * 0.8), zorder=4,
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=1.6))
            ax.text(x * f * 0.72, y * f * 0.72, f'{np.hypot(x, y):.0f}"', color=c, fontsize=6,
                    ha="center", va="center", zorder=4, path_effects=halo)
            continue
        ax.plot(x, y, m, ms=s, mfc="none", mec=c, mew=1.2, zorder=4)
    ax.plot(0, 0, "x", ms=9, mew=1.8, color=SN_COLOR, zorder=5)
    ax.set_xlim(half, -half)  # east left
    ax.set_ylim(-half, half)  # north up
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    fg = "white" if args.dark else "black"
    # 5" scale bar, bottom right; compass, top left
    x1 = -half + 0.08 * size
    ax.plot([x1 + args.scalebar, x1], [-half + 0.08 * size] * 2, color=fg, lw=2)
    ax.text(x1 + args.scalebar / 2, -half + 0.10 * size, f'{args.scalebar:g}"', color=fg, path_effects=halo,
            ha="center", va="bottom", fontsize=7)
    cx, cy, L = half - 0.22 * size, half - 0.22 * size, 0.12 * size  # arrow tips must stay inside
    ax.annotate("", (cx, cy + L), (cx, cy), arrowprops=dict(arrowstyle="-|>", color=fg, lw=1))
    ax.annotate("", (cx + L, cy), (cx, cy), arrowprops=dict(arrowstyle="-|>", color=fg, lw=1))
    ax.text(cx, cy + L * 1.15, "N", color=fg, path_effects=halo, ha="center", va="bottom", fontsize=6)
    ax.text(cx + L * 1.15, cy, "E", color=fg, path_effects=halo, ha="right", va="center", fontsize=6)
    ax.set_title(_title(case), fontsize=7)


def _title(case) -> str:
    title = str(case["case_id"])
    label = case.get("label")
    if isinstance(label, str) and label.strip():
        title += "\n" + textwrap.fill(textwrap.shorten(label, 96, placeholder="..."), 38)
    return title


def render(cases: pd.DataFrame, images: dict, errors: dict, args):
    names = marker_names(cases.columns)
    styles, n_extra = {}, 0
    for n in names:
        styles[n] = marker_style(n, n_extra)
        n_extra += n not in MARKER_STYLES
    present = [n for n in names if cases[f"{n}_ra"].notna().any()]
    handles = [Line2D([], [], ls="none", marker="x", ms=7, mew=1.8, color=SN_COLOR, label="SN")]
    handles += [Line2D([], [], ls="none", marker=styles[n][0], ms=min(styles[n][2], 10), mfc="none",
                       mec=styles[n][1], mew=1.4, label=styles[n][3]) for n in present]

    per_page = args.per_page or len(cases)
    pages = [cases.iloc[i:i + per_page] for i in range(0, len(cases), per_page)]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pdf = PdfPages(out) if out.suffix.lower() == ".pdf" else None
    for ip, page in enumerate(pages):
        ncol = min(args.ncol, len(page))
        nrow = int(np.ceil(len(page) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.75 * nrow + 0.45), squeeze=False,
                                 layout="constrained")
        for ax, (_, case) in zip(axes.flat, page.iterrows()):
            key = case["case_id"]
            draw_panel(ax, case, images.get(key), errors.get(key), names, styles, args.size, args)
        for ax in axes.flat[len(page):]:
            ax.axis("off")
        fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), fontsize=7, frameon=False)
        fig.suptitle(f"{TITLE} {args.band}-band, {args.size:g}\" cutouts"
                     + (f", ellipses at {args.ellipse_scale:g} x the drawn axes" if args.ellipse_scale != 1 else "")
                     + (f" (page {ip + 1}/{len(pages)})" if len(pages) > 1 else ""), fontsize=9)
        if pdf:
            pdf.savefig(fig, dpi=args.dpi)
        else:
            dest = out if len(pages) == 1 else out.with_name(f"{out.stem}_p{ip + 1}{out.suffix}")
            fig.savefig(dest, dpi=args.dpi)
        plt.close(fig)
    if pdf:
        pdf.close()


# ----------------------------------------------------------------------- CLI

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cases", help="csv with case_id, sn_ra, sn_dec, [label], and <name>_ra/<name>_dec pairs")
    p.add_argument("--output", "-o", default="gallery.pdf", help="output .pdf (multi-page) or .png")
    p.add_argument("--band", default="r", choices=list("ugrizy"))
    p.add_argument("--size", type=float, default=30.0, help="cutout side in arcsec")
    p.add_argument("--ncol", type=int, default=4)
    p.add_argument("--per-page", type=int, default=24, help="panels per page (0 = all on one page)")
    p.add_argument("--cache-dir", default="cutout_cache", help="FITS + ObsCore cache directory")
    p.add_argument("--stretch", default="asinh", choices=["asinh", "sqrt", "linear"])
    p.add_argument("--pmin", type=float, default=20.0, help="lower clip percentile (20 ~ just below sky)")
    p.add_argument("--pmax", type=float, default=99.7, help="upper clip percentile")
    p.add_argument("--dark", action="store_true",
                   help="light sources on a dark sky (default is inverted greyscale, where black/grey markers read)")
    p.add_argument("--scalebar", type=float, default=5.0, help="scale bar length in arcsec")
    p.add_argument("--workers", type=int, default=4, help="parallel cutout requests")
    p.add_argument("--min-interval", type=float, default=0.5, help="seconds between request starts")
    p.add_argument("--limit", type=int, default=0, help="render only the first N cases")
    p.add_argument("--sample", type=int, default=0, help="render N random cases (see --seed)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--lookup", help="comma-separated marker names given as <name>_objectId; "
                   "positions are looked up in the DP0.2 Object table")
    p.add_argument("--ellipse-scale", type=float, default=1.0, help="draw <name>_a/_b ellipses scaled by this")
    p.add_argument("--tier", default="dp02", choices=sorted(TIERS), help="which release's coadds")
    p.add_argument("--token-file", help="token file with both read:tap and read:image scopes")
    p.add_argument("--tap-token-file")
    p.add_argument("--image-token-file")
    args = p.parse_args(argv)
    global OBSCORE_TABLE, OBS_WHERE, COADD_SUBTYPE, TITLE
    OBSCORE_TABLE, OBS_WHERE, COADD_SUBTYPE, TITLE = TIERS[args.tier]

    lookup = [n for n in (args.lookup or "").split(",") if n]
    cases = pd.read_csv(args.cases, dtype={f"{n}_objectId": "Int64" for n in lookup})
    missing = {"case_id", "sn_ra", "sn_dec"} - set(cases.columns)
    if missing:
        raise SystemExit(f"{args.cases} lacks columns {sorted(missing)}")
    if args.sample and args.sample < len(cases):
        cases = cases.sample(args.sample, random_state=args.seed)
    if args.limit:
        cases = cases.head(args.limit)
    cases = cases.reset_index(drop=True)
    if cases["case_id"].duplicated().any():
        raise SystemExit("case_id values must be unique")

    cache_dir = Path(args.cache_dir)
    tap_token, image_token = load_tokens(args.token_file, args.tap_token_file, args.image_token_file)
    if lookup:
        cases = lookup_positions(cases, lookup, tap_token)
    images, errors = {}, {}
    try:
        patches = patch_table(args.band, tap_token, cache_dir)
    except Exception as e:  # noqa: BLE001 - render empty panels with the reason
        patches, errors = None, {k: f"ObsCore lookup: {e}" for k in cases["case_id"]}
    if patches is not None:
        session = _session(image_token)
        limiter = RateLimiter(args.min_interval)

        def job(case):
            t0 = time.monotonic()
            content = fetch_cutout(float(case.sn_ra), float(case.sn_dec), args.band, args.size,
                                   patches, session, cache_dir, limiter)
            return read_image(content), time.monotonic() - t0

        t_all = time.monotonic()
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {case.case_id: pool.submit(job, case) for case in cases.itertuples()}
            for key, fut in futures.items():
                try:
                    images[key], dt = fut.result()
                    print(f"  {key}: {dt:.1f} s")
                except Exception as e:  # noqa: BLE001 - one bad case must not kill the gallery
                    errors[key] = f"{type(e).__name__}: {e}"
                    print(f"  {key}: FAILED {errors[key]}")
        print(f"{len(images)}/{len(cases)} cutouts in {time.monotonic() - t_all:.1f} s")
    render(cases, images, errors, args)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
