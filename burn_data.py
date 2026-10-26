"""
Sentinel-2 images of the rice belt north-west of Lahore, one per day with a pass, on one 20 m grid.

Area: about 60 x 60 km over Sheikhupura, Nankana Sahib and Muridke (73.58-74.22 E, 31.35-31.89 N),
in UTM zone 43N, corners on whole kilometres. 20 m is the pixel of the shortwave infrared bands that
burn indices need; the 10 m bands are averaged to it.

Bands kept: blue, green, red (for true colour), narrow NIR (B8A, 20 m), SWIR-1 (B11) and SWIR-2 (B12),
plus a valid-pixel mask from the scene classification (SCL). Valid: vegetation, bare soil, water,
unclassified and also "dark area" pixels (class 2), because fresh char is exactly that; not valid: no
data, saturated, cloud shadow, clouds and thin cirrus. Since processing baseline 04.00 (2022) L2A
values carry an offset of +1000, removed here.

    python burn_data.py 2026-10-01 2026-10-27       # download and cache every day with a pass
"""
import os
import sys
from collections import defaultdict
from datetime import date as Date

import numpy as np
import planetary_computer
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds
from rasterio.windows import Window

CRS = "EPSG:32643"
LONLAT = (73.58, 31.35, 74.22, 31.89)
RES = 20
BANDS = ["B02", "B03", "B04", "B8A", "B11", "B12"]
VALID_SCL = (2, 4, 5, 6, 7)
CACHE = "cache"
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def grid():
    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform
    w, s, e, n = LONLAT
    xs, ys = zip(*(to_utm(x, y) for x in (w, e) for y in (s, n)))
    x0, x1 = np.floor(min(xs) / 1000) * 1000, np.ceil(max(xs) / 1000) * 1000
    y0, top = np.floor(min(ys) / 1000) * 1000, np.ceil(max(ys) / 1000) * 1000
    return from_origin(x0, top, RES, RES), (int((top - y0) // RES), int((x1 - x0) // RES))


def catalog():
    return pystac_client.Client.open(STAC, modifier=planetary_computer.sign_inplace)


def pass_days(start, end):
    """{day: [items]} for every day with a Sentinel-2 pass over the area."""
    days = defaultdict(list)
    for i in catalog().search(collections=["sentinel-2-l2a"], bbox=LONLAT, datetime=f"{start}/{end}").items():
        days[str(i.datetime.date())].append(i)
    return dict(sorted(days.items()))


def _read(href, tr, shape, resampling):
    dst = np.zeros(shape, dtype="float32")
    with rasterio.open(href) as src:
        b = transform_bounds(CRS, src.crs, *rasterio.transform.array_bounds(*shape, tr))
        win = rasterio.windows.from_bounds(b[0] - 100, b[1] - 100, b[2] + 100, b[3] + 100, src.transform)
        try:
            win = win.intersection(Window(0, 0, src.width, src.height)).round_offsets().round_lengths()
        except rasterio.errors.WindowError:
            return dst
        reproject(src.read(1, window=win).astype("float32"), dst, src_transform=src.window_transform(win),
                  src_crs=src.crs, dst_transform=tr, dst_crs=CRS, resampling=resampling, src_nodata=0, dst_nodata=0)
    return dst


def day(d, items=None):
    """(6, H, W) uint16 reflectance x 10000 and the valid mask for one day, cached in cache/."""
    path = os.path.join(CACHE, f"s2_{d}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return z["refl"], z["valid"]
    items = items if items is not None else pass_days(d, d)[d]
    tr, shape = grid()
    refl = np.zeros((len(BANDS),) + shape, dtype="uint16")
    scl = np.zeros(shape, dtype="uint8")
    for item in sorted(items, key=lambda i: i.properties["eo:cloud_cover"]):
        empty = scl == 0
        if not empty.any():
            break
        s = _read(item.assets["SCL"].href, tr, shape, Resampling.nearest)
        take = empty & (s > 0)
        if not take.any():
            continue
        offset = 1000 if item.properties.get("s2:processing_baseline", "00.00") >= "04.00" else 0
        for k, band in enumerate(BANDS):
            res = Resampling.average if band in ("B02", "B03", "B04") else Resampling.nearest
            v = _read(item.assets[band].href, tr, shape, res)
            refl[k][take] = np.clip(v[take] - offset, 1, 60000).astype("uint16")
        scl[take] = s[take]
    valid = np.isin(scl, VALID_SCL)
    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(path, refl=refl, valid=valid)
    return refl, valid


if __name__ == "__main__":
    start, end = sys.argv[1], sys.argv[2]
    for d, items in pass_days(start, end).items():
        refl, valid = day(d, items)
        print(f"{d}: {len(items)} scenes, {valid.mean():.0%} of the area valid", flush=True)
