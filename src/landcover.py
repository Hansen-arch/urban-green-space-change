"""Land-cover attribution for the Pekanbaru vegetation change analysis.

Uses ESA WorldCover 2021 (10 m) to say what kind of land lost or gained
vegetation. The product is a static map from the middle of the study period, so
a label describes the land around 2021 rather than its 2019 state; this is
stated wherever the attribution is reported.
"""

from __future__ import annotations

import ee

import config

# WorldCover values in dict order, and the compact 1..N codes we reduce on.
WORLDCOVER_VALUES = list(config.LANDCOVER_CLASSES)
INDEX = {value: i + 1 for i, value in enumerate(WORLDCOVER_VALUES)}
LABELS = {i + 1: name for i, name in enumerate(config.LANDCOVER_CLASSES.values())}
LABELS[0] = "other"  # snow, mangrove, moss: nominally absent here, kept explicit


def get_landcover(aoi: ee.Geometry) -> ee.Image:
    """WorldCover remapped to compact class codes (0 = other)."""
    raw = ee.Image(config.LANDCOVER_ASSET).select(config.LANDCOVER_BAND).clip(aoi)
    return (
        raw.remap(WORLDCOVER_VALUES, list(INDEX.values()), 0)
        .rename("lc_code")
        .clip(aoi)
    )


def cross_image(change_class: ee.Image, lc_code: ee.Image, mask: ee.Image) -> ee.Image:
    """Encode change class x land cover as one integer: class * 10 + lc_code.

    A single grouped reduction then yields the whole cross-tab, instead of one
    reduceRegion per class pair.
    """
    return change_class.multiply(10).add(lc_code).updateMask(mask).rename("cross")


def parse_groups(groups: list[dict]) -> dict[int, dict[int, float]]:
    """Turn [{'code': 21, 'sum': 1234.5}] (m2) into {change_class: {lc: ha}}."""
    table: dict[int, dict[int, float]] = {1: {}, 2: {}, 3: {}}
    for group in groups or []:
        code = int(group["code"])
        change_class, lc_code = divmod(code, 10)
        if change_class in table:
            table[change_class][lc_code] = group["sum"] / 1e4
    return table
