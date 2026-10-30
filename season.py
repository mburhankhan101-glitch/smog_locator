"""
A whole burning season (1 October to 31 December), year against year, by tehsil.

Detection is step 1's, run one Sentinel-2 day at a time. A season is about 40 days of about 90 MB
each, too much to keep, so days are read without caching and the state of the detection (each
pixel's last clear look, and what was flagged) is saved after every day in cache/season_<year>.npz.
A run that stops picks up where it left off. Days of step 1 already in cache/ are reused.

    python season.py run 2025              # detect through the 2025 season (resumable)
    python season.py run 2026              # the 2026 season so far; run again as it goes on
    python season.py run 2026 --fresh      # start the season over (after a fix to the detection)
    python season.py report                # tables and the figure, both years
"""
import json
import os
import sys
import time
from datetime import date as Date, timedelta

import numpy as np
import pandas as pd

import burn_data
import burns

STATE = os.path.join(burn_data.CACHE, "season_{}.npz")
DETECTION_BANDS = ("B03", "B8A", "B11", "B12")     # green (water mask), NIR, SWIR-1, SWIR-2


def season_days(year):
    end = min(Date(year, 12, 31), Date.today())
    return list(burn_data.pass_days(f"{year}-10-01", str(end)))


def load(year):
    """(state, rows, days) saved so far, or a fresh start."""
    path = STATE.format(year)
    if not os.path.exists(path):
        return burns.start(), [], []
    z = np.load(path, allow_pickle=False)
    s = {k: z[k] for k in burns.start()}
    return s, json.loads(str(z["rows"])), list(z["days"])


def save(year, s, rows, days):
    path = STATE.format(year)
    np.savez_compressed(path + ".tmp.npz", rows=json.dumps(rows), days=np.array(days), **s)
    os.replace(path + ".tmp.npz", path)


def run(year, budget_min=None, fresh=False):
    t0 = time.time()
    days = season_days(year)
    s, rows, done = (burns.start(), [], []) if fresh else load(year)
    if done != days[:len(done)]:
        sys.exit(f"saved state has days {done[:3]}..., not a start of {days[:3]}...; move {STATE.format(year)} away")
    print(f"{year}: {len(days)} Sentinel-2 days, {len(done)} done", flush=True)
    for k in range(len(done), len(days)):
        if budget_min and time.time() - t0 > budget_min * 60:
            print("budget used up; run again to continue", flush=True)
            return
        d = days[k]
        for attempt in range(3):
            try:
                refl, valid = burn_data.day(d, cache=False, bands=DETECTION_BANDS)
                break
            except Exception as e:                      # a dropped connection: search and sign again
                print(f"{d}: {type(e).__name__}, retrying", flush=True)
                time.sleep(20)
        else:
            raise RuntimeError(f"{d}: could not be read")
        rows.append(burns.step(s, days, k, refl, valid))
        save(year, s, rows, days[:k + 1])
        r = rows[-1]
        print(f"{d}: valid {r['valid_share']:.0%}, compared {r['compared_share']:.0%}, burned {r['newly_burned_ha']:,.0f} ha, "
              f"possibly {r['newly_possibly_burned_ha']:,.0f} ha  ({(time.time() - t0) / 60:.0f} min)", flush=True)
    print("season done so far", flush=True)


def tehsils():
    """Tehsil index of every pixel (-1: none) and the names, from data/tehsils.geojson (geoBoundaries)."""
    from pyproj import Transformer
    from rasterio.features import rasterize
    from shapely.geometry import shape
    from shapely.ops import transform
    tr, grid_shape = burn_data.grid()
    to_utm = Transformer.from_crs("EPSG:4326", burn_data.CRS, always_xy=True).transform
    features = json.load(open(os.path.join("data", "tehsils.geojson"), encoding="utf-8"))["features"]
    names = [f["properties"]["tehsil"] for f in features]
    shapes = [(transform(to_utm, shape(f["geometry"])), k) for k, f in enumerate(features)]
    return rasterize(shapes, out_shape=grid_shape, transform=tr, fill=-1, dtype="int16"), names


def year_data(year):
    """Everything the report needs for one season."""
    s, rows, days = load(year)
    r = np.load(os.path.join(burns.OUT, f"rice_{year}.npz"))
    info = pd.read_csv(os.path.join(burns.OUT, f"rice_{year}_groups.csv"))
    table = pd.DataFrame(rows)
    table["date"] = pd.to_datetime(table.day)
    return {"year": year, "days": days, "table": table, "burned_on": s["burned_on"], "possible_on": s["possible_on"],
            "group": r["group"], "rice": r["rice"], "info": info}


def until(y, month_day):
    """Burned pixels first seen on or before month-day (MM-DD) of the year."""
    last = max([k for k, d in enumerate(y["days"]) if d[5:] <= month_day], default=-1)
    return (y["burned_on"] >= 0) & (y["burned_on"] <= last)


def report():
    """Year against year on cropland (ESA WorldCover), which does not depend on the radar: in 2025
    Sentinel-1 orbit 100 covered only 68% of the area, so the rice maps are compared only where both
    years have one (`common`). Rice groups are compared within each year."""
    import rice
    ys = [year_data(year) for year in (2025, 2026) if os.path.exists(STATE.format(year))]
    latest = min(y["days"][-1][5:] for y in ys)                 # compare both years up to the same date
    teh, names = tehsils()
    crop = rice.worldcover() == rice.CROPLAND
    common = np.logical_and.reduce([y["group"] >= 0 for y in ys])
    ha = lambda m: m.sum() * burns.HA
    # 1. same dates, both years (and each year's whole run)
    rows = []
    for y in ys:
        for up_to in sorted({latest, y["days"][-1][5:]}):
            b = until(y, up_to)
            t = y["table"][y["table"].day.str[5:] <= up_to]
            rice_c = y["rice"] & common
            rows.append({"season": y["year"], "up_to": up_to, "sentinel2_days": len(t),
                         "mean_compared_share": t.compared_share.mean(), "burned_ha": ha(b),
                         "burned_on_cropland_ha": ha(b & crop), "share_of_cropland_burned": (b & crop).sum() / crop.sum(),
                         "rice_ha_common_area": ha(rice_c),
                         "share_of_rice_burned_common_area": (b & rice_c).sum() / rice_c.sum()})
    same = pd.DataFrame(rows)
    # 2. by tehsil (those at least 50 km2 inside the area): cropland burned, same dates
    rows = []
    for k, name in enumerate(names):
        m = teh == k
        if ha(m) < 5000:
            continue
        row = {"tehsil": name, "inside_km2": ha(m) / 100, "cropland_ha": ha(m & crop)}
        for y in ys:
            b = until(y, latest) & m & crop
            row[f"burned_on_cropland_ha_{y['year']}"] = ha(b)
            row[f"share_of_cropland_burned_{y['year']}"] = b.sum() / max((m & crop).sum(), 1)
        for y in ys:
            if y["days"][-1][5:] > latest:
                row[f"burned_on_cropland_ha_{y['year']}_whole_run"] = ha((y["burned_on"] >= 0) & m & crop)
        rows.append(row)
    by_tehsil = pd.DataFrame(rows).sort_values("inside_km2", ascending=False)
    # 3. by rice group: when each planting group burned (within a year, each year's own rice map)
    rows = []
    for y in ys:
        doy = np.array([(Date.fromisoformat(d) - Date(y["year"], 10, 1)).days for d in y["days"]])
        for k, g in y["info"].iterrows():
            if not g.rice:
                continue
            m = y["group"] == k
            b = m & (y["burned_on"] >= 0)
            when = doy[y["burned_on"][b]]
            rows.append({"season": y["year"], "group": g.group, "flooded": g.darkest, "area_ha": ha(m),
                         f"burned_by_{latest}": (until(y, latest) & m).sum() / m.sum(),
                         "burned_whole_run": b.sum() / m.sum(),
                         "median_burn_seen": f"{Date(y['year'], 10, 1) + timedelta(days=int(np.median(when))):%m-%d}"
                                             if len(when) else ""})
    by_group = pd.DataFrame(rows)
    out = os.path.join(burns.OUT, "season")
    same.round(4).to_csv(out + "_same_dates.csv", index=False)
    by_tehsil.round(4).to_csv(out + "_tehsils.csv", index=False)
    by_group.round(4).to_csv(out + "_rice_groups.csv", index=False)
    for y in ys:
        t = y["table"]
        on_crop = np.bincount(y["burned_on"][crop & (y["burned_on"] >= 0)], minlength=len(y["days"])) * burns.HA
        t["newly_burned_on_cropland_ha"] = on_crop
        t["burned_on_cropland_so_far_ha"] = np.cumsum(on_crop)
        t.drop(columns="date").round(4).to_csv(os.path.join(burns.OUT, f"season_{y['year']}.csv"), index=False)
    figure(ys, latest, same, by_tehsil, by_group, out + ".png")
    for name, t in (("same dates", same), ("tehsils", by_tehsil), ("rice groups", by_group)):
        print(f"\n{name}\n" + t.round(3).to_string(index=False))


def figure(ys, latest, same, by_tehsil, by_group, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colour = {2025: "#7f8c8d", 2026: "#c0392b"}
    fig = plt.figure(figsize=(20, 13), layout="constrained")
    g = fig.add_gridspec(2, 3)
    ax = [fig.add_subplot(g[0, :2]), fig.add_subplot(g[0, 2]), fig.add_subplot(g[1, 0]),
          fig.add_subplot(g[1, 1]), fig.add_subplot(g[1, 2])]
    for y in ys:
        t = y["table"]
        x = [Date(2000, d.month, d.day) for d in t.date]           # one calendar for both years
        ax[0].plot(x, t.burned_on_cropland_so_far_ha, color=colour[y["year"]], marker="o", lw=2, label=str(y["year"]))
        ax[0].scatter(x, [0] * len(x), marker="|", s=80, color=colour[y["year"]],
                      alpha=np.clip(t.compared_share.values, 0.15, 1))
    ax[0].axvline(Date(2000, int(latest[:2]), int(latest[3:])), color="black", lw=0.8, ls=":")
    ax[0].set_ylabel("cropland seen burned so far (ha)")
    ax[0].legend(fontsize=10)
    ax[0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    ax[0].set_title("Cropland burned so far, through the season (each dot is a Sentinel-2 day; ticks at the bottom are "
                    "fainter when clouds or fog hid most of the area)\ndotted line: the last date both seasons have",
                    fontsize=10)
    w = 0.38
    for i, y in enumerate(ys):
        t = y["table"]
        week = ((t.date - pd.Timestamp(y["year"], 10, 1)).dt.days // 7).values
        wk = pd.Series(t.newly_burned_on_cropland_ha.values).groupby(week).sum()
        ax[1].bar(wk.index + (i - 0.5) * w, wk.values, w, color=colour[y["year"]], label=str(y["year"]))
    ax[1].set_xticks(range(13), [f"{(pd.Timestamp(2000, 10, 1) + pd.Timedelta(days=7 * k)):%d %b}" for k in range(13)],
                     rotation=60, fontsize=8)
    ax[1].set_ylabel("cropland newly seen burned (ha)")
    ax[1].legend(fontsize=9)
    ax[1].set_title("Per week (week starting)", fontsize=10)
    tt = by_tehsil.head(8)
    yy = np.arange(len(tt))
    for i, y in enumerate(ys):
        ax[2].barh(yy + (i - 0.5) * w, tt[f"share_of_cropland_burned_{y['year']}"] * 100, w, color=colour[y["year"]],
                   label=str(y["year"]))
    ax[2].set_yticks(yy, [f"{n} ({km:,.0f} km²)" for n, km in zip(tt.tehsil, tt.inside_km2)], fontsize=9)
    ax[2].invert_yaxis()
    ax[2].set_xlabel(f"share of cropland burned, 1 Oct to {latest[3:]}/{latest[:2]} (%)")
    ax[2].legend(fontsize=9)
    ax[2].set_title("By tehsil, same dates (km²: part inside the study area)", fontsize=10)
    for y in ys:
        gr = by_group[by_group.season == y["year"]]
        x = [Date(2000, int(f[5:7]), int(f[8:10])) for f in gr.flooded]
        ax[3].plot(x, gr[f"burned_by_{latest}"] * 100, color=colour[y["year"]], marker="o", lw=2,
                   label=f"{y['year']}, by {latest[3:]}/{latest[:2]}")
        if y["days"][-1][5:] > latest:
            ax[3].plot(x, gr.burned_whole_run * 100, color=colour[y["year"]], marker="s", lw=1, ls="--",
                       label=f"{y['year']}, whole run to {y['days'][-1][8:]}/{y['days'][-1][5:7]}")
        for xi, v, a in zip(x, gr[f"burned_by_{latest}"] * 100, gr.area_ha):
            ax[3].annotate(f"{a / 1000:.0f}k ha", (xi, v), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax[3].xaxis.set_major_locator(matplotlib.dates.DayLocator(bymonthday=[1, 15]))
    ax[3].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    ax[3].set_xlabel("when the group's paddies were flooded (from radar)")
    ax[3].set_ylabel("share of the group burned (%)")
    ax[3].legend(fontsize=9)
    ax[3].set_title("Planted earlier, burned earlier? Rice groups from step 3, each year", fontsize=10)
    ax[4].axis("off")
    lines = ["Same dates, 1 October to " + f"{latest[3:]}/{latest[:2]}", ""]
    for r in same[same.up_to == latest].itertuples():
        lines.append(f"{r.season}: {r.burned_on_cropland_ha:,.0f} ha of cropland burned ({r.share_of_cropland_burned:.1%})")
        lines.append(f"      {r.sentinel2_days} Sentinel-2 days, {r.mean_compared_share:.0%} of the area compared on average")
    for r in same[same.up_to != latest].itertuples():
        lines += ["", f"{r.season}, whole run to {r.up_to[3:]}/{r.up_to[:2]}: {r.burned_on_cropland_ha:,.0f} ha of cropland "
                      f"({r.share_of_cropland_burned:.1%})"]
    lines += ["", "Cropland: ESA WorldCover 2021, the same map for both years.",
              "The rice maps differ: in 2025 the radar orbit covered 68% of the area."]
    ax[4].text(0, 0.95, "\n".join(lines), va="top", fontsize=11)
    fig.suptitle("The burning season, year against year: rice belt north-west of Lahore", fontsize=13)
    fig.savefig(out, dpi=80)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        budget = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else None
        run(int(sys.argv[2]), budget, "--fresh" in sys.argv)
    elif sys.argv[1] == "report":
        report()
