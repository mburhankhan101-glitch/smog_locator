"""
A weekly burn-scar map for the public: farmland burned per 1 km cell and per tehsil, 2026 against 2025,
and the numbers behind "why fire points undercount". Writes runs/map/ (data, tehsil table, the page).

With the environment variable SMOG_AREA=districts, after the weekly season run and `districts.py fires`:

    python burn_map.py

The map is deliberately coarse. Punjab's enforcement teams act on satellite fire alerts, and these scars
are not checked on the ground, so the public page shows hectares per 1 km cell (about 100 ha of land,
many fields) and per tehsil, never single fields.

Weeks run from 1 October in steps of 7 days, as in the README. A scar is counted in the week of the first
Sentinel-2 image that shows it; it may have burned a few days earlier, between two clear looks.

Tehsils are geoBoundaries ADM3 (data/tehsils.geojson). Only tehsils with at least MIN_INSIDE of their area
inside the Sentinel-2 grid are reported; for those partly outside, the numbers cover the part inside.
"""
import json
import os
import sys
from datetime import date as Date, timedelta

import numpy as np
import pandas as pd

import burn_data
import burns
import districts
import fires
import rice
import season

CELL_PX = 50                       # 50 x 20 m = 1 km
MIN_INSIDE = 0.7                   # report a tehsil if at least this share of it is inside the grid
SIMPLIFY_M = 120                   # outline simplification for the page
DISTRICT = {"Sheikhupura": "Sheikhupura", "Ferozewala": "Sheikhupura", "Muridke": "Sheikhupura",
            "Safdarabad": "Sheikhupura", "Sharak Pur": "Sheikhupura", "Nankana Sahab": "Nankana Sahib",
            "Shah Kot": "Nankana Sahib", "Shangla Hill": "Nankana Sahib", "Lahore City": "Lahore",
            "Lahore Cantt": "Lahore", "Kot Radha Kishan": "Kasur", "Pattoki": "Kasur", "Kasur": "Kasur",
            "Kamoke": "Gujranwala", "Nowshera Virkhan": "Gujranwala", "Gujranwala": "Gujranwala",
            "Jaranwala": "Faisalabad", "Pindi Bhattian": "Hafizabad", "Hafizabad": "Hafizabad"}
OUT = os.path.join(burns.OUT, "map")
TEMPLATE = "burn_map.html"


def week_of(day, year):
    return (Date.fromisoformat(day) - Date(year, 10, 1)).days // 7


def outlines(names):
    """Each tehsil clipped to the grid and simplified, as SVG path data in km from the grid's top left,
    with a label point and the share of the tehsil inside the grid."""
    from pyproj import Transformer
    from shapely.geometry import box, shape
    from shapely.ops import transform, unary_union
    tr, (rows, cols) = burn_data.grid()
    x0, top = tr.c, tr.f
    frame = box(x0, top - rows * burn_data.RES, x0 + cols * burn_data.RES, top)
    to_utm = Transformer.from_crs("EPSG:4326", burn_data.CRS, always_xy=True).transform
    features = json.load(open(os.path.join("data", "tehsils.geojson"), encoding="utf-8"))["features"]
    polys = {f["properties"]["tehsil"]: transform(to_utm, shape(f["geometry"])) for f in features}

    def path(geom):
        geom = geom.simplify(SIMPLIFY_M)
        parts = getattr(geom, "geoms", [geom])
        d = []
        for p in parts:
            for ring in [p.exterior, *p.interiors]:
                pts = [f"{(x - x0) / 1000:.2f} {(top - y) / 1000:.2f}" for x, y in ring.coords]
                d.append("M" + "L".join(pts) + "Z")
        return "".join(d)

    out = {}
    for name in names:
        p = polys[name]
        inside = p.intersection(frame)
        label = inside.representative_point() if not inside.is_empty else None
        out[name] = {"inside": inside.area / p.area, "path": path(inside) if not inside.is_empty else "",
                     "label": [round((label.x - x0) / 1000, 2), round((top - label.y) / 1000, 2)] if label else None}
    by_district = {}
    for name, d in DISTRICT.items():
        if name in polys:
            by_district.setdefault(d, []).append(polys[name])
    district_paths = {d: path(unary_union(ps).intersection(frame)) for d, ps in by_district.items()
                      if not unary_union(ps).intersection(frame).is_empty}
    return out, district_paths, [cols * burn_data.RES / 1000, rows * burn_data.RES / 1000]


def cells(y, crop, weeks):
    """Burned farmland per 1 km cell and week (first image showing the scar): [[cell, week, ha x 10], ...]."""
    b = y["burned_on"]
    rows, cols = b.shape
    week_of_image = np.array([week_of(d, y["year"]) for d in y["days"]])
    m = (b >= 0) & crop
    rr, cc = np.nonzero(m)
    cell = (rr // CELL_PX) * (cols // CELL_PX) + cc // CELL_PX
    wk = week_of_image[b[rr, cc]]
    keep = wk < weeks
    key = cell[keep].astype(np.int64) * weeks + wk[keep]
    k, n = np.unique(key, return_counts=True)
    return [[int(c), int(w), int(round(v * burns.HA * 10))] for c, w, v in zip(k // weeks, k % weeks, n)]


def tehsil_rows(ys, tehsil, names, crop, info, latest):
    """Per reported tehsil: farmland, burned (2026 to date, 2025 to the same date and whole season), fire
    points, and the share of burned farmland VIIRS missed."""
    y26, y25 = ys[2026], ys[2025]
    k26, k25 = len(y26["days"]) - 1, districts.last_index(y25, latest[5:])
    rows = []
    for t, name in enumerate(names):
        if info[name]["inside"] < MIN_INSIDE:
            continue
        m = (tehsil == t) & crop
        r = {"tehsil": name, "district": DISTRICT.get(name, ""), "inside_share": round(info[name]["inside"], 3),
             "farmland_ha": m.sum() * burns.HA}
        for key, y, k in (("2026", y26, k26), ("2025_same_date", y25, k25), ("2025_season", y25, len(y25["days"]) - 1)):
            b = (y["burned_on"] >= 0) & (y["burned_on"] <= k) & m
            compared = y["compared"] & b
            up_to = y["days"][k]
            f = y["fires"]
            r[f"burned_ha_{key}"] = b.sum() * burns.HA
            r[f"share_{key}"] = b.sum() / max(m.sum(), 1)
            r[f"fire_points_{key}"] = int(((f.district == t) & (f.day <= up_to)).sum())
            r[f"missed_{key}"] = 1 - (y["hit"] & compared).sum() / compared.sum() if compared.sum() > 50 else np.nan
        wk = np.array([week_of(d, 2026) for d in y26["days"]])
        per_image = np.bincount(y26["burned_on"][m & (y26["burned_on"] >= 0)], minlength=len(wk)) * burns.HA
        r["weekly_2026"] = [round(float(per_image[wk == w].sum()), 1) for w in range(wk.max() + 1)]
        rows.append(r)
    return pd.DataFrame(rows)


def area_totals(ys, mask, crop, latest):
    """The same numbers over all reported tehsils together, plus weekly series for both seasons."""
    out, y26, y25 = {}, ys[2026], ys[2025]
    m = mask & crop
    for key, y, k in (("2026", y26, len(y26["days"]) - 1), ("2025_same_date", y25, districts.last_index(y25, latest[5:])),
                      ("2025_season", y25, len(y25["days"]) - 1)):
        b = (y["burned_on"] >= 0) & (y["burned_on"] <= k) & m
        compared = y["compared"] & b
        f = y["fires"]
        out[key] = {"up_to": y["days"][k], "burned_ha": b.sum() * burns.HA, "share": b.sum() / m.sum(),
                    "fire_points": int(((f.district >= 0) & (f.day <= y["days"][k])).sum()),
                    "missed": 1 - (y["hit"] & compared).sum() / compared.sum()}
    out["farmland_ha"] = m.sum() * burns.HA
    weekly = {}
    for year, y in ys.items():
        nweeks = week_of(y["days"][-1], year) + 1
        wk = np.array([week_of(d, year) for d in y["days"]])
        per_image = np.bincount(y["burned_on"][m & (y["burned_on"] >= 0)], minlength=len(y["days"])) * burns.HA
        fp = y["fires"][y["fires"].district >= 0]
        fp_week = np.bincount([week_of(d, year) for d in fp.day], minlength=nweeks)[:nweeks]
        weekly[str(year)] = [{"week": w, "start": str(Date(year, 10, 1) + timedelta(days=7 * w)),
                              "burned_ha": round(float(per_image[wk == w].sum()), 1), "fire_points": int(fp_week[w])}
                             for w in range(nweeks)]
    out["weekly"] = weekly
    return out


def main():
    if burn_data.AREA != "districts":
        sys.exit("set SMOG_AREA=districts first")
    tehsil, names = season.tehsils()
    info, district_paths, size_km = outlines(names)
    crop = rice.worldcover() == rice.CROPLAND
    reported = [t for t, n in enumerate(names) if info[n]["inside"] >= MIN_INSIDE]
    mask = np.isin(tehsil, reported)
    ys = {year: districts.one_season(year, np.where(mask, tehsil, -1)) for year in districts.SEASONS}
    latest = ys[2026]["days"][-1]
    table = tehsil_rows(ys, tehsil, names, crop, info, latest)
    totals = area_totals(ys, mask, crop, latest)
    sizes = districts.patches(list(ys.values()), np.where(mask, 0, -1), crop)
    weeks = {year: week_of(y["days"][-1], year) + 1 for year, y in ys.items()}
    images = [{"day": r["day"], "clear": round(r["valid_share"], 3)} for r in ys[2026]["table"].to_dict("records")]
    data = {
        "updated": str(Date.today()), "latest_image": latest, "size_km": size_km, "cell_km": CELL_PX * burn_data.RES / 1000,
        "grid": {"cols": ys[2026]["burned_on"].shape[1] // CELL_PX, "rows": ys[2026]["burned_on"].shape[0] // CELL_PX},
        "weeks": {str(y): n for y, n in weeks.items()},
        "cells": {str(y): cells(ys[y], crop & mask, weeks[y]) for y in ys},      # only inside reported tehsils
        "tehsils": [{"name": n, "district": DISTRICT.get(n, ""), "reported": info[n]["inside"] >= MIN_INSIDE,
                     "inside": round(info[n]["inside"], 3), "path": info[n]["path"], "label": info[n]["label"]}
                    for n in names if info[n]["path"]],
        "districts": district_paths, "table": json.loads(table.round(4).to_json(orient="records")),
        "totals": totals, "images": images, "patches": json.loads(sizes.round(4).to_json(orient="records")),
        "viirs_hours": fires.VIIRS_HOURS,
    }
    os.makedirs(OUT, exist_ok=True)
    table.drop(columns="weekly_2026").round(3).to_csv(os.path.join(OUT, "tehsils.csv"), index=False)
    blob = json.dumps(data, separators=(",", ":"))
    open(os.path.join(OUT, "burn_map.json"), "w", encoding="utf-8").write(blob)
    page = open(TEMPLATE, encoding="utf-8").read().replace("/*DATA*/null", blob)
    open(os.path.join(OUT, "index.html"), "w", encoding="utf-8").write(page)
    t = totals
    print(f"reported tehsils: {len(reported)}; farmland {t['farmland_ha']:,.0f} ha")
    for key in ("2026", "2025_same_date", "2025_season"):
        print(f"{key:15s} to {t[key]['up_to']}: burned {t[key]['burned_ha']:,.0f} ha ({t[key]['share']:.1%}), "
              f"fire points {t[key]['fire_points']:,}, VIIRS missed {t[key]['missed']:.0%}")
    print(table[["tehsil", "district", "inside_share", "farmland_ha", "burned_ha_2026", "share_2026",
                 "burned_ha_2025_same_date", "fire_points_2026", "missed_2026"]].round(3).to_string(index=False))
    print(f"page: {os.path.join(OUT, 'index.html')} ({len(page) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
