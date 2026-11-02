"""
Smoke over each district from Sentinel-5P, once a day at about 13:30 local (the time of the VIIRS
fire pass).

    UV aerosol index   (354/388 nm) rises with particles that absorb ultraviolet light, which smoke
                       does and most haze and dust much less; it also works over bright ground and
                       above thin cloud. About 0 is clean or non-absorbing, above 1 a clear smoke layer.
    carbon monoxide    total column (mmol/m2): a gas that burning gives off, which lives for weeks,
                       so it builds up over the whole region rather than marking single fires.

Step 5 first used aerosol optical depth from Sentinel-3 (aod.py), but it was missing on 25 of 28
days in October 2026 and its retrievals drop the haziest scenes. Sentinel-5P near-real-time granules
are used; pixels below the recommended quality (QA) are left out. A district's value is the mean of
the good pixels whose centre lies inside it, kept when at least MIN_SHARE of its pixels are good.

    python smoke.py 2025-10-01 2025-12-31     # one small CSV per granule in cache/s5p/
"""
import io
import os
import sys
import time

import numpy as np
import pandas as pd
import planetary_computer
import requests
import shapely

import aod
import fires

COLLECTION = "sentinel-5p-l2-netcdf"
PRODUCTS = {"aer-ai": ("aerosol_index_354_388", 0.8, 1, "aerosol_index"),           # variable, QA, scale, name
            "co": ("carbonmonoxide_total_column", 0.5, 1000, "co")}
MIN_SHARE = 0.3
FOLDER = os.path.join("cache", "s5p")


def one_granule(item, polygons, tries=3):
    """Rows (one per district) for one granule."""
    import xarray as xr
    variable, qa, scale, name = PRODUCTS[item.properties["s5p:product_name"]]
    href = next(iter(item.assets.values())).href.split("?")[0]
    for k in range(tries):
        try:
            r = requests.get(planetary_computer.sign(href), timeout=900)       # signed now: links signed at search expire
            r.raise_for_status()
            break
        except requests.RequestException:
            if k == tries - 1:
                raise
            time.sleep(20)
    ds = xr.open_dataset(io.BytesIO(r.content), engine="h5netcdf", group="PRODUCT")
    lat, lon = ds.latitude.values[0], ds.longitude.values[0]
    value, good_qa = ds[variable].values[0] * scale, ds.qa_value.values[0] >= qa
    w, s, e, n = shapely.union_all(list(polygons.values())).bounds
    near = (lon > w) & (lon < e) & (lat > s) & (lat < n)
    rows = []
    for district, poly in polygons.items():
        inside = near.copy()
        inside[near] = shapely.contains_xy(poly, lon[near], lat[near])
        ok = inside & good_qa & np.isfinite(value)
        enough = inside.sum() > 0 and ok.sum() >= MIN_SHARE * inside.sum()
        rows.append({"granule": item.id, "time": item.datetime.isoformat(), "district": district, "product": name,
                     "pixels": int(inside.sum()), "good": int(ok.sum()),
                     "value": float(value[ok].mean()) if enough else np.nan})
    return rows


def fetch(start, end):
    os.makedirs(FOLDER, exist_ok=True)
    polygons = aod.districts()
    bbox = shapely.union_all(list(polygons.values())).bounds
    items = [i for i in fires.s3_items(start, end, collection=COLLECTION, bbox=bbox)
             if i.properties.get("s5p:product_name") in PRODUCTS and i.properties.get("s5p:processing_mode") == "NRTI"]
    print(f"{len(items)} Sentinel-5P granules, {start} to {end}", flush=True)
    for k, item in enumerate(items):
        out = os.path.join(FOLDER, item.id + ".csv")
        if os.path.exists(out):
            continue
        rows = pd.DataFrame(one_granule(item, polygons))
        rows.to_csv(out, index=False)
        print(f"{item.datetime:%Y-%m-%d %H:%M} {rows['product'].iloc[0]}: "
              + ", ".join(f"{r.district} {r.value:.2f} ({r.good}/{r.pixels})" for r in rows.itertuples())
              + f"  [{k + 1}/{len(items)}]", flush=True)


def load():
    """One row per day and district: mean aerosol index and CO of the day's granules (local days)."""
    import glob
    f = pd.concat([pd.read_csv(p) for p in glob.glob(os.path.join(FOLDER, "*.csv"))], ignore_index=True)
    f["day"] = (pd.to_datetime(f.time, utc=True, format="ISO8601") + fires.LOCAL).dt.strftime("%Y-%m-%d")
    return f.pivot_table(index=["day", "district"], columns="product", values="value", aggfunc="mean").reset_index()


if __name__ == "__main__":
    fetch(sys.argv[1], sys.argv[2])
