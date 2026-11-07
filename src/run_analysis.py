"""Pekanbaru vegetation change 2019 vs 2025.

Builds one cloud-masked median NDVI composite per period, difference them into
a change image, then quantifies everything with server-side reduceRegion so the
numbers in the README come from Earth Engine rather than from the downloaded
preview rasters.

Outputs
-------
outputs/data/stats.json              all headline numbers
outputs/data/area_summary.csv        loss / stable / gain areas
outputs/data/address_aoi.geojson     the boundary used, for the figures
outputs/data/ndvi_composites.tif     NDVI for both periods
outputs/data/ndvi_change.tif         dNDVI plus the analysis mask
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import ee

import config
import landcover
import s2_utils
from s2_utils import classify_change, init_ee

MAX_PIXELS = int(1e10)
TILE_SCALE = 4


def reduce(image: ee.Image, geometry: ee.Geometry, reducer: ee.Reducer) -> dict:
    """reduceRegion on the project grid, with generous limits for a whole city."""
    return image.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=config.SCALE,
        crs=config.CRS,
        maxPixels=MAX_PIXELS,
        tileScale=TILE_SCALE,
        bestEffort=False,
    ).getInfo()


def require(values: dict, keys: list[str]) -> dict:
    """Fail loudly instead of writing a null into stats.json."""
    missing = [k for k in keys if values.get(k) is None]
    if missing:
        raise RuntimeError(
            f"reduceRegion returned no value for {missing}. That usually means the "
            "band was fully masked over the AOI - check the SCL mask and dates."
        )
    return values


def annual_context(aoi: ee.Geometry) -> dict:
    """Mean NDVI for every year between the two periods.

    Supports one claim in the interpretation: whether 2019 is a typical year or
    an unusually low one. Same SCL mask and median composite as the headline
    analysis, at CONTEXT_SCALE, without the minimum-observation rule.
    """
    years: dict[str, dict] = {}
    for year in range(config.EARLY_YEAR, config.RECENT_YEAR + 1):
        collection = s2_utils.get_collection(aoi, f"{year}-01-01", f"{year + 1}-01-01")
        ndvi = collection.select("ndvi").median().clip(aoi)
        mean = ndvi.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=aoi,
            scale=config.CONTEXT_SCALE,
            crs=config.CRS,
            maxPixels=MAX_PIXELS,
            tileScale=TILE_SCALE,
        ).getInfo()
        years[str(year)] = {
            "n_scenes_used": collection.size().getInfo(),
            "mean_ndvi": mean.get("ndvi"),
        }
        print(f"  context {year}: mean NDVI {years[str(year)]['mean_ndvi']:.3f}")
    return {
        "scale_m": config.CONTEXT_SCALE,
        "crs": config.CRS,
        "cloud_mask": f"SCL in {config.SCL_VALID}",
        "min_observations": 0,
        "note": (
            "Yearly median-composite mean NDVI. Same SCL mask as the headline "
            "analysis but no minimum-observation or permanent-water filter, so these "
            "values are a context check and are not comparable one-for-one with the "
            "period means in stats.json."
        ),
        "years": years,
    }


def download(image: ee.Image, region: ee.Geometry, path: Path) -> Path:
    """Fetch a small GeoTIFF straight from Earth Engine (no Drive round trip)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    url = image.getDownloadURL(
        {
            "region": region,
            "scale": config.SCALE,
            "crs": config.CRS,
            "format": "GEO_TIFF",
            "filePerBand": False,
        }
    )
    urllib.request.urlretrieve(url, path)
    return path


def main() -> int:
    init_ee()
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    config.FIG_DIR.mkdir(parents=True, exist_ok=True)

    feature, aoi = s2_utils.get_aoi()
    city_area_ha = feature.geometry().area(1).divide(1e4).getInfo()
    print(f"AOI     : {config.GAUL_CITY_VALUE} ({config.GAUL_ASSET})")
    print(f"Area    : {city_area_ha:,.0f} ha | grid {config.SCALE} m {config.CRS}")

    # ---------------------------------------------------------------- composites
    composites: dict[str, dict] = {}
    for key, period in config.PERIODS.items():
        raw = (
            ee.ImageCollection(config.S2_COLLECTION)
            .filterBounds(aoi)
            .filterDate(period["start"], period["end"])
        )
        collection = s2_utils.get_collection(aoi, period["start"], period["end"])
        n_raw = raw.size().getInfo()
        n_used = collection.size().getInfo()
        print(f"{period['label']}  : {n_raw} scenes over AOI, {n_used} after "
              f"CLOUDY_PIXEL_PERCENTAGE < {config.MAX_CLOUDY_PIXEL_PCT}")

        ndvi, valid, water = s2_utils.build_composite(collection, aoi)
        composites[key] = {
            "period": period,
            "ndvi": ndvi,
            "valid": valid,
            "water": water,
            "n_raw": n_raw,
            "n_used": n_used,
        }

    early, recent = composites["early"], composites["recent"]

    # ------------------------------------------------------- change and masks
    mask = s2_utils.analyse_mask(
        early["valid"], recent["valid"], early["water"], recent["water"], aoi
    )
    change = recent["ndvi"].subtract(early["ndvi"]).rename("ndvi_change").clip(aoi)
    change_valid = change.updateMask(mask)
    classes = classify_change(change, mask)

    # -------------------------------------------------------------- statistics
    stats_img = ee.Image.cat(
        [
            early["ndvi"].updateMask(mask).rename("ndvi_early"),
            recent["ndvi"].updateMask(mask).rename("ndvi_recent"),
            change_valid.rename("ndvi_change"),
        ]
    )
    reducer = (
        ee.Reducer.mean()
        .combine(ee.Reducer.stdDev(), sharedInputs=True)
        .combine(ee.Reducer.median(), sharedInputs=True)
        .combine(ee.Reducer.minMax(), sharedInputs=True)
        .combine(ee.Reducer.count(), sharedInputs=True)
    )
    s = require(
        reduce(stats_img, aoi, reducer),
        ["ndvi_early_mean", "ndvi_recent_mean", "ndvi_change_mean"],
    )
    print(f"Mean NDVI: {s['ndvi_early_mean']:.3f} -> {s['ndvi_recent_mean']:.3f} "
          f"({s['ndvi_change_mean']:+.3f})")

    # Areas per class, summed from ee.Image.pixelArea() for an exact hectare count.
    px = ee.Image.pixelArea()
    area_img = ee.Image.cat(
        [
            px.updateMask(classes.eq(1)).rename("loss"),
            px.updateMask(classes.eq(2)).rename("stable"),
            px.updateMask(classes.eq(3)).rename("gain"),
        ]
    ).divide(1e4)
    a = reduce(area_img, aoi, ee.Reducer.sum())
    areas_ha = {
        "loss": a.get("loss") or 0.0,
        "stable": a.get("stable") or 0.0,
        "gain": a.get("gain") or 0.0,
    }
    areas_ha["total"] = sum(areas_ha.values())
    print(f"Areas   : loss {areas_ha['loss']:,.0f} ha | stable {areas_ha['stable']:,.0f} ha "
          f"| gain {areas_ha['gain']:,.0f} ha")

    # How much do the areas depend on the +/- threshold?
    sensitivity = []
    for t in config.SENSITIVITY_THRESHOLDS:
        # Band names cannot contain '.', so 0.05 -> 0p05 in the reducer output.
        tag = str(t).replace(".", "p")
        bands = ee.Image.cat(
            [
                px.updateMask(change_valid.lt(-t)).rename(f"loss_{tag}"),
                px.updateMask(change_valid.lte(t).And(change_valid.gte(-t))).rename(f"stable_{tag}"),
                px.updateMask(change_valid.gt(t)).rename(f"gain_{tag}"),
            ]
        ).divide(1e4)
        r = reduce(bands, aoi, ee.Reducer.sum())
        sensitivity.append(
            {
                "threshold": t,
                "loss_ha": r.get(f"loss_{tag}") or 0.0,
                "stable_ha": r.get(f"stable_{tag}") or 0.0,
                "gain_ha": r.get(f"gain_{tag}") or 0.0,
            }
        )

    # ---------------------------------------------------- land-cover attribution
    lc_code = landcover.get_landcover(aoi)
    cross = landcover.cross_image(classes, lc_code, mask)
    grouped = reduce(
        ee.Image.cat([ee.Image.pixelArea().rename("area"), cross]),
        aoi,
        ee.Reducer.sum().group(groupField=1, groupName="code"),
    ).get("groups")
    by_change = landcover.parse_groups(grouped)

    attribution: dict[str, dict] = {}
    for code, label in sorted(landcover.LABELS.items()):
        parts = {cls: by_change[cls].get(code, 0.0) for cls in (1, 2, 3)}
        if not any(parts.values()):
            continue
        attribution[label] = {
            "loss_ha": parts[1],
            "stable_ha": parts[2],
            "gain_ha": parts[3],
            "net_ha": parts[3] - parts[1],
        }
    attribution = dict(
        sorted(attribution.items(), key=lambda kv: kv[1]["loss_ha"] + kv[1]["gain_ha"], reverse=True)
    )
    print("Attribution by 2021 land cover (loss / gain, ha):")
    for label, row in list(attribution.items())[:4]:
        print(f"   {label:12s} {row['loss_ha']:8,.0f} / {row['gain_ha']:8,.0f}")

    # Observation density, so cloud gaps are visible instead of implied.
    obs_img = ee.Image.cat(
        [
            early["valid"].rename("obs_early"),
            recent["valid"].rename("obs_recent"),
        ]
    )
    obs = reduce(obs_img.updateMask(mask), aoi, ee.Reducer.mean().combine(
        ee.Reducer.min(), sharedInputs=True))

    # ------------------------------------------- year-by-year context series
    print("Annual context series")
    context = annual_context(aoi)
    (config.DATA_DIR / "annual_ndvi_context.json").write_text(json.dumps(context, indent=2))

    # ------------------------------------------------------------------ export
    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"name": config.GAUL_CITY_VALUE, "source": config.GAUL_ASSET},
                "geometry": aoi.getInfo(),
            }
        ],
    }
    (config.DATA_DIR / "aoi_pekanbaru.geojson").write_text(json.dumps(geojson))

    early_label, recent_label = early["period"]["label"], recent["period"]["label"]
    download(
        ee.Image.cat(
            [
                early["ndvi"].updateMask(mask).rename(f"ndvi_{early_label}"),
                recent["ndvi"].updateMask(mask).rename(f"ndvi_{recent_label}"),
            ]
        ),
        aoi,
        config.DATA_DIR / "ndvi_composites.tif",
    )
    download(
        ee.Image.cat(
            [
                change_valid.rename("ndvi_change"),
                mask.toByte().rename("analysis_mask"),
            ]
        ),
        aoi,
        config.DATA_DIR / "ndvi_change.tif",
    )
    download(
        lc_code.toByte().rename("landcover"),
        aoi,
        config.DATA_DIR / "landcover_worldcover.tif",
    )

    # ------------------------------------------------------------------- write
    payload = {
        "meta": {
            "city": config.GAUL_CITY_VALUE,
            "gaul_asset": config.GAUL_ASSET,
            "aoi_field": f"{config.GAUL_CITY_FIELD} = {config.GAUL_CITY_VALUE}",
            "collection": config.S2_COLLECTION,
            "crs": config.CRS,
            "scale_m": config.SCALE,
            "cloud_mask": f"SCL in {config.SCL_VALID}",
            "max_cloudy_pixel_pct": config.MAX_CLOUDY_PIXEL_PCT,
            "change_threshold": config.CHANGE_THRESHOLD,
            "water_fraction": config.WATER_FRACTION,
            "min_observations": config.MIN_OBSERVATIONS,
            "gee_project": config.GEE_PROJECT,
        },
        "city_area_ha": city_area_ha,
        "periods": {
            "early": {
                **early["period"],
                "n_scenes_aoi": early["n_raw"],
                "n_scenes_used": early["n_used"],
                "mean_ndvi": s["ndvi_early_mean"],
                "median_ndvi": s["ndvi_early_median"],
                "std_ndvi": s["ndvi_early_stdDev"],
                "min_ndvi": s["ndvi_early_min"],
                "max_ndvi": s["ndvi_early_max"],
                "n_pixels": s["ndvi_early_count"],
            },
            "recent": {
                **recent["period"],
                "n_scenes_aoi": recent["n_raw"],
                "n_scenes_used": recent["n_used"],
                "mean_ndvi": s["ndvi_recent_mean"],
                "median_ndvi": s["ndvi_recent_median"],
                "std_ndvi": s["ndvi_recent_stdDev"],
                "min_ndvi": s["ndvi_recent_min"],
                "max_ndvi": s["ndvi_recent_max"],
                "n_pixels": s["ndvi_recent_count"],
            },
        },
        "change": {
            "mean": s["ndvi_change_mean"],
            "median": s["ndvi_change_median"],
            "std": s["ndvi_change_stdDev"],
            "min": s["ndvi_change_min"],
            "max": s["ndvi_change_max"],
            "n_pixels": s["ndvi_change_count"],
        },
        "coverage": {
            "analysed_land_ha": areas_ha["total"],
            "analysed_land_pct": areas_ha["total"] / city_area_ha * 100,
            "mean_valid_obs_early": obs.get("obs_early_mean"),
            "mean_valid_obs_recent": obs.get("obs_recent_mean"),
            "min_valid_obs_early": obs.get("obs_early_min"),
            "min_valid_obs_recent": obs.get("obs_recent_min"),
        },
        "areas_ha": areas_ha,
        "areas_pct": {k: (v / areas_ha["total"] * 100) for k, v in areas_ha.items() if k != "total"},
        "threshold_sensitivity": sensitivity,
        "landcover_attribution": {
            "asset": config.LANDCOVER_ASSET,
            "class_field": f"remapped from {config.LANDCOVER_BAND}",
            "by_class": attribution,
            "note": (
                "WorldCover 2021 is a static map from the middle of the study period, "
                "so classes describe the land around 2021, not its 2019 state. Tree "
                "plantations fall under 'tree cover', so this cannot separate estate "
                "from natural forest."
            ),
        },
        "net_vegetation_area_change_ha": areas_ha["gain"] - areas_ha["loss"],
    }

    (config.DATA_DIR / "stats.json").write_text(json.dumps(payload, indent=2))

    rows = ["class,area_ha,area_pct_of_analysed_land"]
    for k in ("loss", "stable", "gain"):
        rows.append(f"{k},{areas_ha[k]:.2f},{payload['areas_pct'][k]:.2f}")
    (config.DATA_DIR / "area_summary.csv").write_text("\n".join(rows) + "\n")

    lc_rows = ["landcover,loss_ha,stable_ha,gain_ha,net_ha"]
    for label, row in attribution.items():
        lc_rows.append(
            f"{label},{row['loss_ha']:.2f},{row['stable_ha']:.2f},"
            f"{row['gain_ha']:.2f},{row['net_ha']:.2f}"
        )
    (config.DATA_DIR / "landcover_attribution.csv").write_text("\n".join(lc_rows) + "\n")

    print()
    print(f"Analysed land : {areas_ha['total']:,.0f} ha "
          f"({payload['coverage']['analysed_land_pct']:.1f}% of the city)")
    print(f"Net change    : {payload['net_vegetation_area_change_ha']:+,.0f} ha "
          f"(gain - loss at |dNDVI| > {config.CHANGE_THRESHOLD})")
    print(f"Valid obs/px  : {payload['coverage']['mean_valid_obs_early']:.1f} "
          f"({early_label}) vs {payload['coverage']['mean_valid_obs_recent']:.1f} ({recent_label})")
    print(f"\nWrote {config.DATA_DIR}/stats.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
