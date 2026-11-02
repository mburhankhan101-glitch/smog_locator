"""
Aerosol optical depth (AOD) over each district, from Sentinel-3, one value per pass.

AOD at 550 nm says how much the particles in the whole column of air, ground to space, dim sunlight:
about 0.1 is clean air, above 1 thick haze. Fine-mode AOD counts only the small particles (smoke,
traffic, industry), not dust. It is a column: the same smoke gives about the same AOD whether it is
mixed through 2 km of afternoon air or squeezed under a 200 m inversion, so it is a signal of how
much haze is around, not of what people breathe.

Sentinel-3 SYNERGY aerosol product (SY_2_AOD), about 4.5 km pixels, from Sentinel-3A and 3B at about
10:00-11:00 local time. Retrievals flagged cloudy, failed, invalid, outlier or zero are left out.
A district's value is the mean of the good pixels whose centre lies inside it, kept only when at
least MIN_SHARE of its pixels are good (otherwise clouds or fog hid it, and the value is missing).

    python aod.py 2025-10-01 2025-12-31     # one small CSV per pass in cache/aod/
"""
import io
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import requests
import shapely
from shapely.geometry import shape

import fires

COLLECTION = "sentinel-3-synergy-aod-l2-netcdf"
BAD = 4 | 256 | 1024 | 2048 | 4096          # cloudy (nadir), AOD zero, uncertainty failed, invalid, outlier
MIN_SHARE = 0.3
FOLDER = os.path.join("cache", "aod")
VARIABLES = {"AOD_550": "aod", "FM_AOD550": "fine_aod", "AAOD550": "absorbing_aod"}


def districts():
    """{district: polygon} from data/districts.geojson."""
    features = json.load(open(os.path.join("data", "districts.geojson"), encoding="utf-8"))["features"]
    return {f["properties"]["district"]: shape(f["geometry"]) for f in features}


def one_pass(item, polygons, tries=3):
    """Rows (one per district) for one Sentinel-3 pass."""
    import xarray as xr
    for k in range(tries):
        try:
            r = requests.get(item.assets["ntc-aod"].href, timeout=600)
            r.raise_for_status()
            break
        except requests.RequestException:
            if k == tries - 1:
                raise
            time.sleep(20)
    ds = xr.open_dataset(io.BytesIO(r.content), engine="h5netcdf")
    lat, lon = ds.latitude.values, ds.longitude.values
    w, s, e, n = shapely.union_all(list(polygons.values())).bounds
    near = (lon > w) & (lon < e) & (lat > s) & (lat < n)
    good = near & ~(ds.aod_quality_flags.fillna(2048).values.astype("int64") & BAD).astype(bool) \
        & np.isfinite(ds.AOD_550.values)
    rows = []
    for name, poly in polygons.items():
        inside = near.copy()
        inside[near] = shapely.contains_xy(poly, lon[near], lat[near])
        ok = inside & good
        row = {"pass": item.id, "time": item.datetime.isoformat(), "district": name,
               "pixels": int(inside.sum()), "good": int(ok.sum())}
        enough = inside.sum() > 0 and ok.sum() >= MIN_SHARE * inside.sum()
        for v, col in VARIABLES.items():
            row[col] = float(np.nanmean(ds[v].values[ok])) if enough else np.nan
        rows.append(row)
    return rows


def fetch(start, end):
    os.makedirs(FOLDER, exist_ok=True)
    polygons = districts()
    bbox = shapely.union_all(list(polygons.values())).bounds
    items = list(fires.s3_items(start, end, collection=COLLECTION, bbox=bbox))
    print(f"{len(items)} Sentinel-3 passes, {start} to {end}", flush=True)
    for k, item in enumerate(items):
        out = os.path.join(FOLDER, item.id + ".csv")
        if os.path.exists(out):
            continue
        rows = pd.DataFrame(one_pass(item, polygons))
        rows.to_csv(out, index=False)
        print(f"{item.datetime:%Y-%m-%d %H:%M} {item.id[:3]}: "
              + ", ".join(f"{r.district} {r.aod:.2f} ({r.good}/{r.pixels})" for r in rows.itertuples())
              + f"  [{k + 1}/{len(items)}]", flush=True)


def load():
    """Every pass fetched so far, one row per pass and district, with the local day."""
    import glob
    f = pd.concat([pd.read_csv(p) for p in glob.glob(os.path.join(FOLDER, "*.csv"))], ignore_index=True)
    f["time"] = pd.to_datetime(f.time, utc=True)
    f["day"] = (f.time + fires.LOCAL).dt.strftime("%Y-%m-%d")
    return f.sort_values("time").reset_index(drop=True)


if __name__ == "__main__":
    fetch(sys.argv[1], sys.argv[2])
