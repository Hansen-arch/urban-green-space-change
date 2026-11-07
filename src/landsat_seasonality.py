"""Cross-sensor check: how much does the seasonal window change the answer?

The main analysis compares two single calendar years of Sentinel-2, which is an
annual window rather than a seasonally matched one. This script asks how much
that choice matters, using Landsat 8/9 as an independent sensor with a much
longer archive: it computes NDVI change for 2016-2018 against 2022-2024 twice,
once with a dry-season window (June to September) and once with full years.

The point is not to reproduce the Sentinel-2 numbers. Landsat, different epochs
and a different resampling mean these are not the same measurement. The point is
to bound the seasonal effect: if the dry-season and annual answers disagree, then
an annual window is a real source of uncertainty in the headline result.

Writes outputs/data/landsat_seasonality.csv.
"""
from __future__ import annotations

import csv
import sys

import ee

import config
import s2_utils

WINDOWS = {
    "dry_season": (("2016-06-01", "2018-09-30"), ("2022-06-01", "2024-09-30")),
    "full_year": (("2016-01-01", "2018-12-31"), ("2022-01-01", "2024-12-31")),
}

LANDSAT = ["LANDSAT/LC08/C02/T1_L2", "LANDSAT/LC09/C02/T1_L2"]
THRESHOLD = 0.03


def composite(start: str, end: str, aoi: ee.Geometry) -> ee.Image:
    """Median clear-sky NDVI from Landsat 8/9 surface reflectance."""
    coll = ee.ImageCollection(LANDSAT[0]).merge(ee.ImageCollection(LANDSAT[1]))
    coll = coll.filterBounds(aoi).filterDate(start, end).filter(
        ee.Filter.lt("CLOUD_COVER", 40)
    )

    def prep(image: ee.Image) -> ee.Image:
        qa = image.select("QA_PIXEL")
        clear = (
            qa.bitwiseAnd(1 << 1).eq(0)
            .And(qa.bitwiseAnd(1 << 3).eq(0))
            .And(qa.bitwiseAnd(1 << 4).eq(0))
        )
        sr = image.select(["SR_B5", "SR_B4"]).multiply(0.0000275).add(-0.2)
        ndvi = sr.normalizedDifference(["SR_B5", "SR_B4"]).rename("ndvi")
        return ndvi.updateMask(clear)

    return ee.ImageCollection(coll.map(prep)).median()


def measure(aoi: ee.Geometry, base: tuple[str, str], recent: tuple[str, str]) -> dict:
    before = composite(base[0], base[1], aoi)
    after = composite(recent[0], recent[1], aoi)
    change = after.subtract(before).rename("ndvi_change")
    stack = change.addBands(change.lt(-THRESHOLD).rename("loss"))
    stack = stack.addBands(change.gt(THRESHOLD).rename("gain"))

    stats = stack.reduceRegion(
        reducer=ee.Reducer.mean().combine(ee.Reducer.median(), sharedInputs=True),
        geometry=aoi,
        scale=30,
        maxPixels=1e10,
        bestEffort=True,
    ).getInfo()

    def pick(kind: str, band: str) -> float:
        # combine() names outputs "<band>_<reducer>"
        for key, value in stats.items():
            if key == f"{band}_{kind}":
                return value
        raise KeyError(f"{band}_{kind} missing from reducer output: {sorted(stats)}")

    return {
        "net_change": pick("mean", "ndvi_change"),
        "median_change": pick("median", "ndvi_change"),
        "pct_loss": 100 * pick("mean", "loss"),
        "pct_gain": 100 * pick("mean", "gain"),
    }


def main() -> int:
    s2_utils.init_ee()
    _, aoi = s2_utils.get_aoi()

    rows = []
    for name, (base, recent) in WINDOWS.items():
        result = measure(aoi, base, recent)
        result["window"] = name
        rows.append(result)
        print(
            f"{name:11s} net {result['net_change']:+.3f}  "
            f"median {result['median_change']:+.3f}  "
            f"loss {result['pct_loss']:.1f}%  gain {result['pct_gain']:.1f}%"
        )

    out = config.DATA_DIR / "landsat_seasonality.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["window", "net_change", "median_change", "pct_loss", "pct_gain"],
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
