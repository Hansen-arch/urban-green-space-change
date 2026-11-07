"""Project configuration: AOI, periods, band and threshold constants.

Everything that a reader might want to change lives here so the analysis and
figure scripts stay free of magic numbers.
"""

from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "outputs" / "data"
FIG_DIR = ROOT / "outputs" / "figures"

# --------------------------------------------------------------------------
# Earth Engine
# --------------------------------------------------------------------------
GEE_PROJECT = "loyal-aileron-473318-c6"

# FAO GAUL 2025 level 2. The 2015 release this used to point at is deprecated
# and, more importantly, the 2025 release renamed ADM2_NAME -> GAUL2_NAME.
GAUL_ASSET = "FAO/GAUL/2025/level2"
GAUL_CITY_FIELD = "GAUL2_NAME"
GAUL_CITY_VALUE = "Kota Pekanbaru"

# --------------------------------------------------------------------------
# Periods
# --------------------------------------------------------------------------
# Sentinel-2 SR (L2A) begins over Pekanbaru on 2018-12-15, so 2019 is the first
# year with a complete record. Override here if you need a different pair.
EARLY_YEAR = 2019
RECENT_YEAR = 2025

PERIODS = {
    "early": {"label": str(EARLY_YEAR), "start": f"{EARLY_YEAR}-01-01", "end": f"{EARLY_YEAR + 1}-01-01"},
    "recent": {"label": str(RECENT_YEAR), "start": f"{RECENT_YEAR}-01-01", "end": f"{RECENT_YEAR + 1}-01-01"},
}

# --------------------------------------------------------------------------
# Sentinel-2 SR
# --------------------------------------------------------------------------
S2_COLLECTION = "COPERNICUS/S2_SR_HARMONIZED"

# Loose scene-level filter only: the real cloud handling is the per-pixel SCL
# mask below, and a strict threshold here would throw away scenes that are
# clear over the city but cloudy elsewhere in the tile.
MAX_CLOUDY_PIXEL_PCT = 80

# Scene Classification Layer classes kept as valid observations.
#   4 vegetation, 5 not-vegetated, 6 water, 11 snow/ice
# Everything else is masked: 0 no-data, 1 saturated/defective, 2 dark area,
# 3 cloud shadow, 7 unclassified, 8/9 cloud medium/high, 10 cirrus.
SCL_VALID = [4, 5, 6, 11]
SCL_WATER = 6

RED_BAND = "B4"
NIR_BAND = "B8"

# --------------------------------------------------------------------------
# Analysis grid
# --------------------------------------------------------------------------
# 20 m, not 10 m: SCL is natively 20 m, so this is the finest grid at which the
# cloud mask is a real observation rather than an interpolation of one.
SCALE = 20
CRS = "EPSG:32647"  # UTM zone 47N, lon ~101.5E -> zone 47 (matches Riau convention)

# Coarser grid for the year-by-year context series. It is a supporting check on
# how unusual the two chosen years are, not a headline number, so 40 m keeps the
# extra six composites cheap.
CONTEXT_SCALE = 40

# Land pixels where |dNDVI| is at or below this are classed as stable.
CHANGE_THRESHOLD = 0.10
# Reported alongside the headline threshold so readers can see sensitivity.
SENSITIVITY_THRESHOLDS = [0.05, 0.10, 0.15]

# A pixel counts as permanent water when SCL class 6 covers more than this
# fraction of its clear observations in *both* periods.
WATER_FRACTION = 0.50

# --------------------------------------------------------------------------
# Land-cover attribution
# --------------------------------------------------------------------------
# ESA WorldCover 2021: complete coverage over the city (99.7%) and a single
# authoritative date. It sits mid-period, which matters for interpreting
# attribution: a pixel cleared in 2020 is labelled by its post-clearing state.
LANDCOVER_ASSET = "ESA/WorldCover/v200/2021"
LANDCOVER_BAND = "Map"

# WorldCover value -> label. The dict order sets the compact code used in the
# cross-tab (position + 1); anything unmapped falls to 0 = "other".
LANDCOVER_CLASSES = {
    10: "tree cover",
    20: "shrubland",
    30: "grassland",
    40: "cropland",
    50: "built-up",
    60: "bare",
    80: "water",
    90: "wetland",
}

# Minimum clear observations per pixel per period. A single clear scene is a
# bad basis for a median (one hazy pass can swing NDVI by 0.5), so those pixels
# are dropped. This costs ~0.7% of the city in the early period.
MIN_OBSERVATIONS = 5
