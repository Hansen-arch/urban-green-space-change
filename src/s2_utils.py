"""Sentinel-2 SR helpers: AOI, SCL cloud masking, median composites, NDVI.

All functions return Earth Engine objects (server-side handles); nothing here
triggers computation. Statistics and downloads live in run_analysis.py.
"""

from __future__ import annotations

import warnings

import ee

import config


def init_ee() -> None:
    """Initialise Earth Engine against the project configured in config.py."""
    warnings.filterwarnings("ignore", category=FutureWarning)
    try:
        ee.Initialize(project=config.GEE_PROJECT)
    except Exception as exc:  # noqa: BLE001 - re-raised with actionable text
        raise RuntimeError(
            "Earth Engine is not authenticated on this machine. Run "
            "`earthengine authenticate` once, then re-run this script."
        ) from exc


def get_aoi() -> tuple[ee.Feature, ee.Geometry]:
    """Return the Kota Pekanbaru feature and its geometry."""
    fc = ee.FeatureCollection(config.GAUL_ASSET).filter(
        ee.Filter.eq(config.GAUL_CITY_FIELD, config.GAUL_CITY_VALUE)
    )
    feature = ee.Feature(fc.first())
    return feature, feature.geometry()


def _mask_scl(image: ee.Image) -> ee.Image:
    """Drop clouds, cloud shadow, cirrus, saturated and no-data pixels.

    Uses the SCL band shipped with L2A: only SCL_VALID classes survive, so the
    mask is a per-pixel, per-scene observation rather than a scene-level guess.
    """
    valid = image.select("SCL").remap(config.SCL_VALID, [1] * len(config.SCL_VALID), 0)
    return image.updateMask(valid)


def _add_ndvi(image: ee.Image) -> ee.Image:
    """Scale SR reflectance and append an NDVI band.

    S2_SR_HARMONIZED is harmonised to pre-2022 processing baseline, so the
    -1000 BOA offset has already been removed and no offset is applied here.
    NDVI is invariant to the 1e-4 scale factor, but scaling keeps the
    reflectance values honest for anything added later.
    """
    scaled = image.select([config.RED_BAND, config.NIR_BAND]).multiply(0.0001)
    ndvi = scaled.normalizedDifference([config.NIR_BAND, config.RED_BAND]).rename("ndvi")
    return (
        image.select("SCL")
        .addBands(scaled)
        .addBands(ndvi)
        .copyProperties(image, ["system:time_start"])
    )


def get_collection(aoi: ee.Geometry, start: str, end: str) -> ee.ImageCollection:
    """Cloud-masked, NDVI-bearing collection bounded to the AOI and period."""
    return (
        ee.ImageCollection(config.S2_COLLECTION)
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", config.MAX_CLOUDY_PIXEL_PCT))
        .map(_mask_scl)
        .map(_add_ndvi)
    )


def build_composite(collection: ee.ImageCollection, aoi: ee.Geometry):
    """Per-pixel median NDVI, valid-observation count and water fraction.

    Returns (ndvi, valid_count, water_fraction), each clipped to the AOI.
    valid_count is the number of clear observations behind each median value,
    which is how we show that no month silently disappeared.
    """
    ndvi_col = collection.select("ndvi")

    ndvi = ndvi_col.median().rename("ndvi").clip(aoi)

    valid_count = (
        ndvi_col.map(lambda img: img.mask().unmask(0).byte())
        .sum()
        .rename("valid_obs")
        .clip(aoi)
    )

    water_count = (
        collection.select("SCL")
        .map(lambda img: img.eq(config.SCL_WATER).unmask(0))
        .sum()
        .rename("water_obs")
    )
    # Guard against zero valid observations before dividing.
    water_fraction = (
        water_count.divide(valid_count.max(1))
        .rename("water_fraction")
        .clip(aoi)
    )

    return ndvi, valid_count, water_fraction


def analyse_mask(
    valid_early: ee.Image,
    valid_recent: ee.Image,
    water_early: ee.Image,
    water_recent: ee.Image,
    aoi: ee.Geometry,
) -> ee.Image:
    """Pixels with enough clear observations in both periods, not water.

    Change needs a valid observation in both periods, otherwise a pixel that
    was simply cloudy in one year would masquerade as a change. On top of that,
    a median built from one or two clear scenes is fragile, so pixels below
    MIN_OBSERVATIONS are dropped. Permanent water is excluded so the mean is a
    land-only number.
    """
    min_obs = config.MIN_OBSERVATIONS
    observed_both = valid_early.gte(min_obs).And(valid_recent.gte(min_obs))
    permanent_water = water_early.gt(config.WATER_FRACTION).And(
        water_recent.gt(config.WATER_FRACTION)
    )
    return observed_both.And(permanent_water.Not()).rename("analysis_mask").clip(aoi)


def classify_change(change: ee.Image, mask: ee.Image) -> ee.Image:
    """Class band: 1 loss, 2 stable, 3 gain (0 / masked elsewhere).

    Thresholds are symmetric around zero so loss and gain are directly
    comparable; SENSITIVITY_THRESHOLDS in the stats section tests how much the
    headline areas depend on this choice.
    """
    t = config.CHANGE_THRESHOLD
    return (
        ee.Image(0)
        .updateMask(mask)
        .where(change.lt(-t), 1)
        .where(change.abs().lte(t), 2)
        .where(change.gt(t), 3)
        .rename("change_class")
    )
