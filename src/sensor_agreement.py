"""Do two sensors agree on where vegetation changed?

The main analysis uses Sentinel-2 for 2019 and 2025. Landsat 8/9 covers a longer
period, so it can be asked the same question over 2016-2018 against 2022-2024.
The two sensors see different resolution, revisits and atmospheric handling, and
they cover different periods, so their magnitudes are not expected to match. What
should match is the spatial pattern: if both put the loss in the same cells, the
Sentinel-2 result is supported by an independent instrument.

Cells are 2 km squares over the area covered by both rasters. A cell is used only
if it holds at least 100 valid pixels in each, so edge slivers do not drive the
correlation.

Writes outputs/data/sensor_agreement.csv and outputs/figures/fig7_sensor_agreement.png.
"""
from __future__ import annotations

import argparse
import csv

import ee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import requests
from rasterio.features import geometry_mask
from rasterio.io import MemoryFile
from scipy import stats
from shapely.geometry import box

import config
import landsat_seasonality as ls
import s2_utils

CELL_M = 2000
MIN_PIXELS = 100


def fetch_landsat_change(path, aoi) -> None:
    """Download Landsat NDVI change at 30 m for the dry-season window."""
    before = ls.composite(*ls.WINDOWS["dry_season"][0], aoi)
    after = ls.composite(*ls.WINDOWS["dry_season"][1], aoi)
    change = after.subtract(before).rename("ndvi_change")

    url = change.getDownloadURL({
        "scale": 30, "crs": config.CRS, "region": aoi,
        "format": "GEO_TIFF", "filePerBand": False,
    })
    response = requests.get(url, timeout=600)
    response.raise_for_status()
    with MemoryFile(response.content) as mf, mf.open() as src:
        profile = src.profile.copy()
        array = src.read(1)
    profile.update(count=1, dtype="float32", compress="deflate")
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype("float32"), 1)
    print(f"wrote {path}")


def read_raster(path, band=1):
    with rasterio.open(path) as src:
        return src.read(band), src.transform, src.bounds


def cell_pairs(sentinel, smask, s_transform, s_bounds,
               landsat, l_transform, l_bounds):
    x0, y0 = max(s_bounds.left, l_bounds.left), max(s_bounds.bottom, l_bounds.bottom)
    x1, y1 = min(s_bounds.right, l_bounds.right), min(s_bounds.top, l_bounds.top)
    xs = np.arange(np.ceil(x0 / CELL_M) * CELL_M, x1, CELL_M)
    ys = np.arange(np.ceil(y0 / CELL_M) * CELL_M, y1, CELL_M)

    rows = []
    for x in xs:
        for y in ys:
            geom = box(x, y, x + CELL_M, y + CELL_M)
            m_s = geometry_mask([geom], sentinel.shape, s_transform, invert=True)
            m_l = geometry_mask([geom], landsat.shape, l_transform, invert=True)
            v_s = sentinel[m_s & (smask > 0)]
            v_l = landsat[m_l]
            if v_s.size < MIN_PIXELS or v_l.size < MIN_PIXELS:
                continue
            rows.append({
                "x": int(x), "y": int(y),
                "sentinel_change": float(v_s.mean()),
                "landsat_change": float(v_l.mean()),
                "sentinel_pixels": int(v_s.size),
                "landsat_pixels": int(v_l.size),
            })
    return rows


def figure(rows, pearson, spearman, path):
    s = np.array([r["sentinel_change"] for r in rows])
    l = np.array([r["landsat_change"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.4, 6))
    ax.scatter(l, s, s=22, alpha=0.75, edgecolor="white", linewidth=0.4, color="#2b6cb0")
    lim = [min(l.min(), s.min()), max(l.max(), s.max())]
    ax.plot(lim, lim, color="black", linewidth=0.8, linestyle="--", label="1:1")
    slope, intercept = np.polyfit(l, s, 1)
    grid = np.linspace(l.min(), l.max(), 50)
    ax.plot(grid, intercept + slope * grid, color="#c53030", linewidth=1.2,
            label=f"fit: {slope:.2f}x + {intercept:.3f}")
    ax.set_xlabel("Landsat NDVI change, 2016-2018 to 2022-2024")
    ax.set_ylabel("Sentinel-2 NDVI change, 2019 to 2025")
    ax.set_title(f"Sensor agreement per 2 km cell (n = {len(rows)})\n"
                 f"Pearson r = {pearson:.2f}, Spearman rho = {spearman:.2f}", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sentinel", default=str(config.DATA_DIR / "ndvi_change.tif"))
    ap.add_argument("--out-csv", default=str(config.DATA_DIR / "sensor_agreement.csv"))
    ap.add_argument("--out-fig", default=str(config.FIG_DIR / "fig7_sensor_agreement.png"))
    args = ap.parse_args()

    landsat_path = config.DATA_DIR / "landsat_change.tif"
    if not landsat_path.exists():
        s2_utils.init_ee()
        _, aoi = s2_utils.get_aoi()
        fetch_landsat_change(str(landsat_path), aoi)

    sentinel, s_transform, s_bounds = read_raster(args.sentinel, 1)
    smask, _, _ = read_raster(args.sentinel, 2)
    landsat, l_transform, l_bounds = read_raster(str(landsat_path), 1)

    rows = cell_pairs(sentinel, smask, s_transform, s_bounds,
                      landsat, l_transform, l_bounds)
    if len(rows) < 10:
        print(f"only {len(rows)} comparable cells, not enough to correlate")
        return 1

    s = np.array([r["sentinel_change"] for r in rows])
    l = np.array([r["landsat_change"] for r in rows])
    pearson = stats.pearsonr(s, l)
    spearman = stats.spearmanr(s, l)

    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure(rows, pearson.statistic, spearman.statistic, args.out_fig)
    print(f"cells {len(rows)}")
    print(f"pearson r = {pearson.statistic:.3f} (p = {pearson.pvalue:.2g})")
    print(f"spearman rho = {spearman.statistic:.3f} (p = {spearman.pvalue:.2g})")
    print(f"mean Sentinel {s.mean():+.3f}, mean Landsat {l.mean():+.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
