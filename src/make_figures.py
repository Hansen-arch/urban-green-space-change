"""Figures for the Pekanbaru vegetation change project.

Reads the GeoTIFFs and stats.json written by run_analysis.py. The statistics in
the figures come from Earth Engine; the rasters are only used for rendering.

Outputs
-------
outputs/figures/fig1_ndvi_before_after.png
outputs/figures/fig2_ndvi_change_map.png
outputs/figures/fig3_change_graphs.png
"""

from __future__ import annotations

import json
import sys

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

import config

GREY = "#eef0f2"
INK = "#1a1a1a"
LOSS = "#b2182b"
STABLE = "#bdbdbd"
GAIN = "#1b7837"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": "#9a9a9a",
        "axes.linewidth": 0.6,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    }
)


# ---------------------------------------------------------------- loading
def load_stats() -> dict:
    return json.loads((config.DATA_DIR / "stats.json").read_text())


def load_boundary():
    return gpd.read_file(config.DATA_DIR / "aoi_pekanbaru.geojson").to_crs(config.CRS)


def load_raster(path):
    """Return (array, extent, transform) with masked pixels set to NaN."""
    with rasterio.open(path) as src:
        arr = src.read().astype("float32")
        extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
    return arr, extent


# ------------------------------------------------------------- map helpers
def style_map(ax, extent, boundary):
    ax.set_facecolor(GREY)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    boundary.boundary.plot(ax=ax, color="#3a3a3a", linewidth=0.9, zorder=5)


def scale_bar(ax, extent, length_m=5000):
    """Scale bar on a white plate, so the label survives dark land cover."""
    x0 = extent[0] + 0.05 * (extent[1] - extent[0])
    y0 = extent[2] + 0.05 * (extent[3] - extent[2])
    ax.add_patch(
        Rectangle((x0, y0), length_m * 1.04, 620, facecolor="white",
                  edgecolor="none", alpha=0.92, zorder=6)
    )
    bar_y = y0 + 170
    ax.plot([x0 + length_m * 0.02, x0 + length_m * 1.02], [bar_y, bar_y],
            color=INK, linewidth=1.6, zorder=7)
    ax.text(
        x0 + length_m * 0.52,
        bar_y + 120,
        f"{length_m / 1000:g} km",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color=INK,
        zorder=7,
    )


def north_arrow(ax, extent):
    x = extent[1] - 0.06 * (extent[1] - extent[0])
    y = extent[3] - 0.14 * (extent[3] - extent[2])
    h = 0.055 * (extent[3] - extent[2])
    ax.annotate(
        "",
        xy=(x, y + h),
        xytext=(x, y),
        arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.4),
        zorder=7,
    )
    ax.text(x, y + h * 1.15, "N", ha="center", va="bottom", fontsize=8, color=INK, zorder=7)


# ------------------------------------------------------------------ fig 1
def figure_before_after(stats, ndvi, extent, boundary, out):
    early_label = stats["periods"]["early"]["label"]
    recent_label = stats["periods"]["recent"]["label"]

    fig = plt.figure(figsize=(10.6, 6.2))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[1, 1, 0.045], wspace=0.05,
        left=0.015, right=0.915, top=0.845, bottom=0.13,
    )
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    cax = fig.add_subplot(gs[0, 2])
    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad(GREY)
    norm = matplotlib.colors.Normalize(vmin=-0.2, vmax=0.9)

    panels = [(axes[0], ndvi[0], early_label), (axes[1], ndvi[1], recent_label)]
    for ax, arr, label in panels:
        plot = np.where(np.isfinite(arr), arr, np.nan)
        im = ax.imshow(plot, extent=extent, origin="upper", cmap=cmap, norm=norm,
                       interpolation="nearest")
        style_map(ax, extent, boundary)
        scale_bar(ax, extent)
        north_arrow(ax, extent)
        mean = stats["periods"]["early" if label == early_label else "recent"]["mean_ndvi"]
        ax.set_title(f"{label}\nmean NDVI {mean:.3f}", fontsize=10.5, color=INK, pad=7)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("NDVI", fontsize=9)
    cbar.ax.tick_params(labelsize=7.5)

    fig.text(
        0.465, 0.965,
        "Vegetation greenness in Kota Pekanbaru, Riau",
        fontsize=13, color=INK, ha="center", va="top",
    )
    fig.text(
        0.465, 0.912,
        f"Sentinel-2 SR, {config.SCALE} m grid, annual median composites, cloud-masked",
        ha="center", fontsize=8.5, color="#4d4d4d",
    )
    fig.text(
        0.015, 0.055,
        "Source: Copernicus Sentinel-2 L2A (COPERNICUS/S2_SR_HARMONIZED, SCL cloud mask) via Google Earth Engine.\n"
        f"Boundary: FAO GAUL 2025 level 2. Projection EPSG:32647 (UTM 47N). Grey = masked or outside AOI.\n"
        f"Pixels shown only where both periods have ≥ {config.MIN_OBSERVATIONS} clear observations.",
        ha="left", va="top", fontsize=6.8, color="#5a5a5a",
    )
    fig.savefig(out)
    plt.close(fig)


# ------------------------------------------------------------------ fig 2
def figure_change(stats, change, mask, extent, boundary, out):
    early = stats["periods"]["early"]["label"]
    recent = stats["periods"]["recent"]["label"]
    d = np.where(mask == 1, change, np.nan)

    fig, ax = plt.subplots(figsize=(7.8, 7.9))
    fig.subplots_adjust(left=0.02, right=0.985, top=0.86, bottom=0.11)
    cmap = plt.get_cmap("RdYlGn").copy()
    cmap.set_bad(GREY)
    norm = TwoSlopeNorm(vmin=-0.4, vcenter=0.0, vmax=0.4)
    im = ax.imshow(d, extent=extent, origin="upper", cmap=cmap, norm=norm, interpolation="nearest")
    style_map(ax, extent, boundary)
    scale_bar(ax, extent)
    north_arrow(ax, extent)

    cbar = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.02, extend="both")
    cbar.set_label(f"ΔNDVI ({recent} − {early})", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5)

    fig.text(
        0.5, 0.965,
        f"Change in vegetation greenness, {early} → {recent}",
        fontsize=12.5, color=INK, ha="center", va="top",
    )
    fig.text(
        0.5, 0.918,
        f"mean ΔNDVI {stats['change']['mean']:+.3f}    ·    "
        f"{stats['areas_ha']['loss']:,.0f} ha lost, {stats['areas_ha']['gain']:,.0f} ha gained "
        f"at |ΔNDVI| > {config.CHANGE_THRESHOLD}",
        ha="center", va="top", fontsize=8.6, color="#4d4d4d",
    )
    fig.text(
        0.02, 0.055,
        "Red = vegetation decline, green = increase. Grey = masked (cloud in both periods, permanent water,\n"
        f"or < {config.MIN_OBSERVATIONS} clear observations) or outside the city boundary.\n"
        "Source: Copernicus Sentinel-2 L2A via Google Earth Engine. Projection EPSG:32647 (UTM 47N).",
        ha="left", va="top", fontsize=6.8, color="#5a5a5a",
    )
    fig.savefig(out)
    plt.close(fig)


# ------------------------------------------------------------------ fig 3
def figure_graphs(stats, change, mask, out):
    areas = stats["areas_ha"]
    pct = stats["areas_pct"]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(10.6, 4.5), gridspec_kw={"width_ratios": [1, 1.35]}
    )

    # -- panel a: area by change class
    labels = ["Loss", "Stable", "Gain"]
    keys = ["loss", "stable", "gain"]
    values = [areas[k] for k in keys]
    colors = [LOSS, STABLE, GAIN]
    bars = ax1.bar(labels, values, color=colors, width=0.62, edgecolor="none")
    for bar, v, k in zip(bars, values, keys):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            v + max(values) * 0.02,
            f"{v:,.0f} ha\n{pct[k]:.1f}%",
            ha="center", va="bottom", fontsize=8.5, color=INK,
        )
    ax1.set_ylabel("Area (ha)", fontsize=9)
    ax1.set_ylim(0, max(values) * 1.22)
    ax1.set_title("Land area by change class", fontsize=10.5, color=INK)
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.tick_params(labelsize=8.5)
    ax1.text(
        0.02, 0.95,
        f"net {stats['net_vegetation_area_change_ha']:+,.0f} ha",
        transform=ax1.transAxes, ha="left", va="top", fontsize=9, color="#4d4d4d",
    )

    # -- panel b: distribution of pixel-level change
    d = change[mask == 1]
    d = d[np.isfinite(d)]
    bins = np.arange(-0.6, 0.601, 0.01)
    ax2.hist(d, bins=bins, color="#4a6fa5", edgecolor="none")
    for t in (config.CHANGE_THRESHOLD, -config.CHANGE_THRESHOLD):
        ax2.axvline(t, color=LOSS, linestyle="--", linewidth=0.9)
    ax2.axvline(0, color="#666666", linewidth=0.7)
    ax2.axvline(stats["change"]["mean"], color=GAIN, linewidth=1.4)
    ax2.set_xlim(-0.6, 0.6)
    ax2.set_xlabel("ΔNDVI (2025 − 2019)", fontsize=9)
    ax2.set_ylabel("Pixels", fontsize=9)
    ax2.set_title("Distribution of pixel-level NDVI change", fontsize=10.5, color=INK)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.tick_params(labelsize=8.5)
    ax2.text(
        0.02, 0.95,
        f"mean {stats['change']['mean']:+.3f}\nmedian {stats['change']['median']:+.3f}\n"
        f"n = {stats['change']['n_pixels']:,} px",
        transform=ax2.transAxes, ha="left", va="top", fontsize=7.8, color="#4d4d4d",
    )
    ax2.legend(
        handles=[
            plt.Line2D([], [], color=LOSS, linestyle="--", linewidth=0.9,
                       label=f"±{config.CHANGE_THRESHOLD} threshold"),
            plt.Line2D([], [], color=GAIN, linewidth=1.4, label="mean change"),
        ],
        loc="upper right", frameon=False, fontsize=7.5, handlelength=1.8,
    )

    fig.suptitle(
        f"NDVI change in Kota Pekanbaru, {stats['periods']['early']['label']} → "
        f"{stats['periods']['recent']['label']}",
        fontsize=12, color=INK, y=1.0,
    )
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def figure_annual_context(context, out):
    """Small context chart: how typical are the two years being compared?"""
    years = context["years"]
    x = [int(y) for y in years]
    y = [years[str(yr)]["mean_ndvi"] for yr in x]

    fig, ax = plt.subplots(figsize=(7.2, 3.5))
    fig.subplots_adjust(left=0.10, right=0.975, top=0.80, bottom=0.155)

    ax.plot(x, y, color="#4a6fa5", linewidth=1.7, zorder=2)
    ax.scatter(x, y, s=30, color="#4a6fa5", zorder=3)

    # Circle the two years the project actually compares.
    for yr, colour in ((x[0], LOSS), (x[-1], GAIN)):
        ax.scatter(
            [yr], [y[x.index(yr)]], s=110, facecolor="white",
            edgecolor=colour, linewidth=1.8, zorder=4,
        )

    for yr, val in zip(x, y):
        ax.annotate(
            f"{val:.3f}", (yr, val), textcoords="offset points",
            xytext=(0, 9), ha="center", fontsize=7.4, color="#4d4d4d",
        )

    peak = int(np.argmax(y))
    ax.annotate(
        f"peak {x[peak]}", (x[peak], y[peak]), textcoords="offset points",
        xytext=(0, 25), ha="center", fontsize=7.6, color=INK,
    )

    ax.set_xticks(x)
    ax.set_xticklabels([str(v) for v in x], fontsize=8.5)
    ax.set_ylim(min(y) - 0.008, max(y) + 0.028)
    ax.set_ylabel("Mean NDVI", fontsize=9)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8.5)

    ax.legend(
        handles=[
            plt.Line2D([], [], marker="o", linestyle="none", markersize=7,
                       markerfacecolor="white", markeredgecolor="#555555",
                       label="years compared"),
        ],
        loc="lower right", frameon=False, fontsize=7.4, handletextpad=0.4,
    )

    fig.text(
        0.5, 0.965,
        "Annual mean NDVI over Kota Pekanbaru, " + f"{x[0]}–{x[-1]}",
        fontsize=11, color=INK, ha="center", va="top",
    )
    fig.text(
        0.10, 0.075,
        f"Median composites at {context['scale_m']} m, SCL-masked, no minimum-observation filter.\n"
        "Context only: not directly comparable with the period means in the table above.",
        ha="left", va="top", fontsize=6.8, color="#5a5a5a",
    )
    fig.savefig(out)
    plt.close(fig)


def figure_landcover_attribution(stats, out):
    """Diverging bars: what kind of land lost or gained vegetation."""
    table = stats["landcover_attribution"]["by_class"]
    labels = list(table)
    loss = [table[k]["loss_ha"] for k in labels]
    gain = [table[k]["gain_ha"] for k in labels]

    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    fig.subplots_adjust(left=0.155, right=0.975, top=0.80, bottom=0.20)

    y = np.arange(len(labels))
    ax.barh(y, [-v for v in loss], color=LOSS, height=0.66, label="loss")
    ax.barh(y, gain, color=GAIN, height=0.66, label="gain")

    span = max(max(loss), max(gain)) * 1.22
    for yi, l, g in zip(y, loss, gain):
        ax.text(-l, yi, f" {l:,.0f}", ha="right", va="center", fontsize=7.2, color=LOSS)
        ax.text(g, yi, f" {g:,.0f}", ha="left", va="center", fontsize=7.2, color=GAIN)

    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(-span, span)
    ax.set_xticks([-10000, -5000, 0, 5000, 10000])
    ax.set_xticklabels(["10,000", "5,000", "0", "5,000", "10,000"], fontsize=8)
    ax.set_xlabel("Area (ha)   ← lost        gained →", fontsize=8.8)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax.legend(loc="lower right", frameon=False, fontsize=8, ncol=2,
              bbox_to_anchor=(1.0, 0.02))

    fig.text(
        0.5, 0.965,
        "Where vegetation was lost and gained, by land cover",
        fontsize=11, color=INK, ha="center", va="top",
    )
    fig.text(
        0.155, 0.905,
        f"ESA WorldCover 2021 · change at |ΔNDVI| > {config.CHANGE_THRESHOLD} · area in hectares",
        ha="left", va="top", fontsize=8.2, color="#4d4d4d",
    )
    fig.text(
        0.155, 0.10,
        "WorldCover is a single map dated mid-period: a class describes the land around 2021, not its 2019 state.\n"
        "Tree plantations are classed as 'tree cover', so this cannot separate estate from natural forest.",
        ha="left", va="top", fontsize=6.8, color="#5a5a5a",
    )
    fig.savefig(out)
    plt.close(fig)


LC_COLOURS = {
    1: "#006400",  # tree cover
    2: "#FFBB22",  # shrubland
    3: "#FFFF4C",  # grassland
    4: "#F096FF",  # cropland
    5: "#FA0000",  # built-up
    6: "#B4B4B4",  # bare
    7: "#0064C8",  # water
    8: "#0096A0",  # wetland
}


def figure_landcover_map(stats, out):
    """WorldCover classes as a map, so the attribution has a geography."""
    arr, extent = load_raster(config.DATA_DIR / "landcover_worldcover.tif")
    band = arr[0]
    band = np.where(band >= 1, band, np.nan)  # 0 = outside the city
    boundary = load_boundary()

    labels = {i + 1: name for i, name in enumerate(config.LANDCOVER_CLASSES.values())}
    codes = sorted(labels)
    cmap = matplotlib.colors.ListedColormap([LC_COLOURS[c] for c in codes])
    cmap.set_bad(GREY)
    norm = matplotlib.colors.BoundaryNorm(np.arange(0.5, len(codes) + 1.5, 1), len(codes))

    fig, ax = plt.subplots(figsize=(7.8, 8.0))
    fig.subplots_adjust(left=0.02, right=0.985, top=0.86, bottom=0.17)
    ax.imshow(band, extent=extent, origin="upper", cmap=cmap, norm=norm,
              interpolation="nearest")
    style_map(ax, extent, boundary)
    scale_bar(ax, extent)
    north_arrow(ax, extent)

    handles = [
        Rectangle((0, 0), 1, 1, facecolor=LC_COLOURS[c], edgecolor="none",
                  label=labels[c])
        for c in codes
    ]
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.035),
        ncol=4, frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.4,
        handletextpad=0.5,
    )

    fig.text(
        0.5, 0.965,
        "Land cover used for attribution, ESA WorldCover 2021",
        fontsize=12, color=INK, ha="center", va="top",
    )
    fig.text(
        0.5, 0.922,
        f"10 m product, shown on the {config.SCALE} m analysis grid over the Kota Pekanbaru boundary",
        ha="center", fontsize=8.2, color="#4d4d4d",
    )
    fig.text(
        0.02, 0.085,
        "Tree plantations are classed as tree cover, which is why cropland covers only ~145 ha of the city.\n"
        "Grey = outside the city boundary.",
        ha="left", va="top", fontsize=6.8, color="#5a5a5a",
    )
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)
    stats = load_stats()
    boundary = load_boundary()

    ndvi, ndvi_extent = load_raster(config.DATA_DIR / "ndvi_composites.tif")
    change_stack, change_extent = load_raster(config.DATA_DIR / "ndvi_change.tif")
    change, mask = change_stack[0], change_stack[1]

    figure_before_after(
        stats, ndvi, ndvi_extent, boundary, config.FIG_DIR / "fig1_ndvi_before_after.png"
    )
    figure_change(
        stats, change, mask, change_extent, boundary, config.FIG_DIR / "fig2_ndvi_change_map.png"
    )
    figure_graphs(stats, change, mask, config.FIG_DIR / "fig3_change_graphs.png")

    context = json.loads((config.DATA_DIR / "annual_ndvi_context.json").read_text())
    figure_annual_context(context, config.FIG_DIR / "fig4_annual_ndvi_context.png")
    figure_landcover_map(stats, config.FIG_DIR / "fig5_landcover_map.png")
    figure_landcover_attribution(
        stats, config.FIG_DIR / "fig6_landcover_attribution.png"
    )

    for path in sorted(config.FIG_DIR.glob("*.png")):
        print(f"wrote {path.relative_to(config.ROOT)} ({path.stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
