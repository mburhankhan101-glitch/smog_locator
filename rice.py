"""
Where rice grew in kharif 2026, and which burn scars lie on it.

Rice in Punjab is transplanted into flooded paddies in late June and July. To Sentinel-1 radar,
standing water is a mirror: the pulse bounces away and VH backscatter is very low. As the rice
grows, the pulse bounces between the stems and the water and back to the satellite, and VH rises
by several dB through August and September. Sugarcane, orchards and fodder are never flooded and
stay flat; villages stay bright; ponds stay dark.

Sentinel-1 RTC (gamma0, terrain corrected) from one orbit only (ORBIT, ascending, about 18:00
local), so every date sees the fields from the same angle: 11 dates from late June to mid October.
VH is averaged to the 20 m grid of step 1 and over 3 x 3 pixels (60 m) in linear power, against
speckle. On cropland (ESA WorldCover 2021, class 40) the 11-date VH curves are grouped by k-means
(K groups, fitted on a random sample). A group is rice when its average curve was flooded (lowest VH
from mid June to mid August at or below FLOOD_DB) and then rose by at least RISE_DB. Rules on
pixels did not work: neither the lowest VH nor the rise splits into two clear humps, so any
pixel cut-off is arbitrary, while the group curves are clean.

    python rice.py                 # writes runs/rice_2026.npz, rice_2026_groups.csv, rice_2026_burns.csv, .png
"""
import os
from datetime import date as Date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rasterio.warp import Resampling
from scipy import ndimage

import burn_data
import burns
import fires

ORBIT = 100
SEASON = ("2026-06-01", "2026-10-28")
FLOOD_WINDOW = ("2026-06-15", "2026-08-15")
FLOOD_DB, RISE_DB = -19.5, 4.0
K = 6
CROPLAND = 40
WORLDCOVER = {10: "trees", 20: "shrubs", 30: "grass", 40: "cropland", 50: "built-up", 60: "bare",
              80: "water", 90: "wetland"}
BURNS = os.path.join(burns.OUT, "burns_2026-10-01_2026-10-27")
FIRES = os.path.join(burns.OUT, "fires_2026-10-01_2026-10-27_points.csv")
OUT = os.path.join(burns.OUT, "rice_2026")


def s1_days():
    """{day: item} for orbit ORBIT over the season, searched a month at a time (long searches time out)."""
    days = {}
    start, end = Date.fromisoformat(SEASON[0]), Date.fromisoformat(SEASON[1])
    for m in range(start.month, end.month + 1):
        a = max(start, Date(start.year, m, 1))
        b = min(end, Date.fromordinal(Date(start.year, m + 1, 1).toordinal() - 1))
        for i in burn_data.catalog().search(collections=["sentinel-1-rtc"], bbox=burn_data.LONLAT,
                                            datetime=f"{a}/{b}T23:59:59Z", method="GET", limit=100).items():
            if i.properties["sat:relative_orbit"] == ORBIT:
                days[str(i.datetime.date())] = i
    return dict(sorted(days.items()))


def s1_day(d, item=None):
    """VV and VH in dB on the 20 m grid for one date, cached as cache/s1_<day>.npz (float16)."""
    path = os.path.join(burn_data.CACHE, f"s1_{d}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return z["vv"].astype("float32"), z["vh"].astype("float32")
    tr, shape = burn_data.grid()
    out = {}
    for pol in ("vv", "vh"):
        lin = burn_data._read(item.assets[pol].href, tr, shape, Resampling.average)
        with np.errstate(divide="ignore"):
            out[pol] = np.where(lin > 0, 10 * np.log10(lin), np.nan).astype("float16")
    os.makedirs(burn_data.CACHE, exist_ok=True)
    np.savez_compressed(path, **out)
    return out["vv"].astype("float32"), out["vh"].astype("float32")


def worldcover():
    """ESA WorldCover 2021 class on the 20 m grid (the most common 10 m class in each pixel)."""
    path = os.path.join(burn_data.CACHE, "worldcover_2021.npy")
    if os.path.exists(path):
        return np.load(path)
    tr, shape = burn_data.grid()
    lc = np.zeros(shape, dtype="uint8")
    for item in burn_data.catalog().search(collections=["esa-worldcover"], bbox=burn_data.LONLAT, method="GET").items():
        if item.properties.get("esa_worldcover:product_version", "2.0.0").startswith("2"):
            v = burn_data._read(item.assets["map"].href, tr, shape, Resampling.mode)
            lc[(lc == 0) & (v > 0)] = v[(lc == 0) & (v > 0)]
    np.save(path, lc)
    return lc


def smooth_db(db, size=3):
    """Mean over size x size pixels in linear power, back to dB."""
    lin = np.where(np.isnan(db), 0, 10 ** (db / 10))
    n = ndimage.uniform_filter((~np.isnan(db)).astype("float32"), size)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10 * np.log10(ndimage.uniform_filter(lin, size) / n)


def groups(days, vh, lc):
    """k-means groups of the cropland VH curves. Returns the group of every pixel (-1: not cropland
    or no data), the average curves, and for each group whether it is rice and when it was darkest.
    Rice groups come first, ordered by flood date, then the others from dark to bright."""
    from sklearn.cluster import KMeans
    crop = (lc == CROPLAND) & np.isfinite(vh).all(0)
    x = vh[:, crop].T
    sample = np.random.default_rng(0).choice(len(x), min(200_000, len(x)), replace=False)
    km = KMeans(K, n_init=4, random_state=0).fit(x[sample])
    centres = km.cluster_centers_
    t = np.array([Date.fromisoformat(d) for d in days])
    fw = (t >= Date.fromisoformat(FLOOD_WINDOW[0])) & (t <= Date.fromisoformat(FLOOD_WINDOW[1]))
    lowest = np.where(fw, centres, np.inf).argmin(1)
    low = centres[np.arange(K), lowest]
    rise = np.array([centres[k, lowest[k]:].max() - low[k] for k in range(K)])
    rice = (low <= FLOOD_DB) & (rise >= RISE_DB)
    order = sorted(range(K), key=lambda k: (not rice[k], lowest[k] if rice[k] else -centres[k].mean()))
    new = np.empty(K, dtype=int)
    new[order] = np.arange(K)
    group = np.full(lc.shape, -1, dtype="int8")
    group[crop] = new[km.predict(x)]
    info = pd.DataFrame({"rice": rice[order], "darkest": [days[lowest[k]] for k in order],
                         "lowest_db": low[order].round(1), "rise_db": rise[order].round(1)})
    return group, centres[order], info


def name(row):
    if row.rice:
        return f"rice, flooded ~{Date.fromisoformat(row.darkest):%d %b}" + (", weaker" if row.lowest_db > -21 else "")
    return "not flooded, dark" if row.lowest_db < -16 else "not flooded, bright"


def missed(points, burned_on, since, days, mask):
    """VIIRS: burned area (inside mask) compared, caught, and the share missed (as in fires.py)."""
    g = points[(points.sensor == "VIIRS") & ~points.persistent].reset_index(drop=True)
    first, last = fires.coverage((g.time.min(), g.time.max()), days)
    flagged = np.where(mask, burned_on, -1)
    hit, _, lo, hi, _ = fires.compare(g, flagged, since, days)
    inside = (flagged >= 0) & (lo >= first) & (hi <= last)
    return inside.sum() * burns.HA, (hit & inside).sum() * burns.HA, 1 - (hit & inside).sum() / max(inside.sum(), 1)


def main():
    days = s1_days()
    for d, item in days.items():
        s1_day(d, item)
    days = list(days)
    vh = np.stack([smooth_db(s1_day(d)[1]) for d in days])
    lc = worldcover()
    group, centres, info = groups(days, vh, lc)
    rice = np.isin(group, np.flatnonzero(info.rice))
    z = np.load(BURNS + ".npz")
    burned_on, possible_on, since, s2_days = z["burned_on"], z["possible_on"], z["since"], list(z["days"])
    burned, possible = burned_on >= 0, possible_on >= 0
    refl, valid = burn_data.day(s2_days[-1])
    r = refl.astype("float32") / 10000
    with np.errstate(invalid="ignore", divide="ignore"):
        ndvi = (r[3] - r[2]) / (r[3] + r[2])
    # per group: area, burns, and greenness on the last Sentinel-2 day (an independent check)
    info.insert(0, "group", [name(row) for row in info.itertuples()])
    info["share_of_cropland"] = [(group == k).sum() / (group >= 0).sum() for k in range(K)]
    info["area_ha"] = [(group == k).sum() * burns.HA for k in range(K)]
    info["burned_ha"] = [(burned & (group == k)).sum() * burns.HA for k in range(K)]
    info["possibly_burned_ha"] = [(possible & (group == k)).sum() * burns.HA for k in range(K)]
    info["burned_share_of_group"] = info.burned_ha / info.area_ha
    info[f"ndvi_{s2_days[-1]}"] = [np.nanmedian(ndvi[(group == k) & valid]) for k in range(K)]
    # where the burns are: rice, other cropland, outside cropland
    places = {"rice": rice, "cropland, not rice": (group >= 0) & ~rice, "cropland, no radar": (lc == CROPLAND) & (group < 0)}
    places.update({f"not cropland: {v}": lc == k for k, v in WORLDCOVER.items() if k != CROPLAND})
    points = pd.read_csv(FIRES)
    points["time"] = pd.to_datetime(points.time, utc=True)
    rows = []
    for where, m in places.items():
        row = {"where": where, "area_ha": m.sum() * burns.HA, "burned_ha": (burned & m).sum() * burns.HA,
               "possibly_burned_ha": (possible & m).sum() * burns.HA}
        if where in ("rice", "cropland, not rice"):
            row["viirs_compared_ha"], row["viirs_caught_ha"], row["viirs_missed_share"] = missed(points, burned_on, since, s2_days, m)
        rows.append(row)
    table = pd.DataFrame(rows)
    table = table[table.burned_ha + table.possibly_burned_ha >= 1]
    np.savez_compressed(OUT + ".npz", group=group, rice=rice, days=np.array(days), centres=centres)
    info.round(3).to_csv(OUT + "_groups.csv", index=False)
    table.round(3).to_csv(OUT + "_burns.csv", index=False)
    figure(days, centres, info, group, lc, burned, possible, table, refl, OUT + ".png")
    print(info.round(2).to_string(index=False))
    print(table.round(2).to_string(index=False))
    print(f"rice: {rice.sum() * burns.HA:,.0f} ha ({rice.sum() / (group >= 0).sum():.0%} of cropland); "
          f"burned on rice {(burned & rice).sum() / burned.sum():.0%}, possibly burned on rice {(possible & rice).sum() / possible.sum():.0%}")


def figure(days, centres, info, group, lc, burned, possible, table, refl, out, zoom_km=4):
    tr, shape = burn_data.grid()
    w_km, h_km = shape[1] * burn_data.RES / 1000, shape[0] * burn_data.RES / 1000
    colours = ["#1a7f37", "#2ea043", "#56d364", "#a8e6a1", "#d4a017", "#8c6d1f"][:K]
    rice_colours = iter(colours[:int(info.rice.sum())])
    other_colours = iter(["#c9a227", "#7a5c12", "#b08850", "#6b4f1d"])
    palette = [next(rice_colours) if r else next(other_colours) for r in info.rice]
    fig = plt.figure(figsize=(20, 13), layout="constrained")
    g = fig.add_gridspec(2, 3)
    ax = [fig.add_subplot(g[0, 0]), fig.add_subplot(g[0, 1]), fig.add_subplot(g[0, 2]),
          fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1]), fig.add_subplot(g[1, 2])]
    t = [Date.fromisoformat(d) for d in days]
    for k, row in enumerate(info.itertuples()):
        ax[0].plot(t, centres[k], color=palette[k], lw=2.5 if row.rice else 1.5, ls="-" if row.rice else "--",
                   marker="o", ms=3, label=f"{row.group} ({row.share_of_cropland:.0%})")
    ax[0].axhline(FLOOD_DB, color="grey", lw=0.8, ls=":")
    ax[0].text(t[0], FLOOD_DB - 0.4, f"flooded: {FLOOD_DB} dB", fontsize=8, color="grey", va="top")
    ax[0].set_ylabel("VH backscatter (dB)")
    ax[0].legend(fontsize=8, loc="lower right")
    ax[0].tick_params(axis="x", labelrotation=45, labelsize=8)
    ax[0].set_title("Average radar curve of each group of cropland pixels\n"
                    "rice: dark when flooded for transplanting, then bright as it grows", fontsize=10)
    rgb = np.ones(shape + (3,)) * 0.85
    rgb[lc == 50] = (0.45, 0.45, 0.45)
    rgb[lc == 80] = (0.3, 0.5, 0.8)
    for k in range(K):
        rgb[group == k] = matplotlib.colors.to_rgb(palette[k])
    ax[1].imshow(rgb[::3, ::3], extent=(0, w_km, 0, h_km), interpolation="nearest")
    ax[1].set_title(f"Rice (greens) and other cropland (browns), kharif 2026\n"
                    f"rice {info.area_ha[info.rice].sum():,.0f} ha, {info.share_of_cropland[info.rice].sum():.0%} of cropland; "
                    f"grey: built-up, blue: water", fontsize=10)
    labels = [w.replace("not cropland: ", "") for w in table["where"]]
    y = np.arange(len(table))
    ax[2].barh(y, table.burned_ha, color="#c0392b", label="burned")
    ax[2].barh(y, table.possibly_burned_ha, left=table.burned_ha, color="#bbbbbb", label="possibly burned")
    for k, r in enumerate(table.itertuples()):
        ax[2].text(r.burned_ha + r.possibly_burned_ha, k, f" {r.burned_ha:,.0f} + {r.possibly_burned_ha:,.0f} ha",
                   va="center", fontsize=8)
    ax[2].set_yticks(y, labels, fontsize=9)
    ax[2].invert_yaxis()
    ax[2].set_xlim(0, (table.burned_ha + table.possibly_burned_ha).max() * 1.35)
    ax[2].set_xlabel("hectares, 1-26 October 2026")
    ax[2].legend(fontsize=8, loc="lower right")
    ax[2].set_title("Where the burn scars lie", fontsize=10)
    ax[3].bar(range(K), info.burned_share_of_group * 100, color=palette)
    short = [gname.replace("rice, flooded ", "rice\n").replace(", weaker", "\nweaker").replace("not flooded, ", "not flooded\n")
             for gname in info.group]
    ax[3].set_xticks(range(K), short, fontsize=8)
    ax[3].set_ylabel("share of the group's area burned (%)")
    ax[3].set_title("Burned so far in October, by group\nfields planted first were harvested and burned first", fontsize=10)
    # zoom: the window with the most burned area, groups underneath, burns outlined
    n = int(zoom_km * 1000 / burn_data.RES)
    dens = np.cumsum(np.cumsum(np.pad(burned.astype("int32"), ((1, 0), (1, 0))), 0), 1)
    box = dens[n:, n:] - dens[:-n, n:] - dens[n:, :-n] + dens[:-n, :-n]
    r0, c0 = np.unravel_index(np.argmax(box), box.shape)
    zz = (slice(r0, r0 + n), slice(c0, c0 + n))
    ax[1].add_patch(plt.Rectangle((c0 * burn_data.RES / 1000, h_km - (r0 + n) * burn_data.RES / 1000), zoom_km, zoom_km,
                                  fill=False, edgecolor="red", lw=1.5))
    ax[4].imshow(rgb[zz], extent=(0, zoom_km, 0, zoom_km), interpolation="nearest")
    ax[4].contour(np.flipud(burned[zz]), [0.5], colors="red", linewidths=0.8, extent=(0, zoom_km, 0, zoom_km))
    ax[4].contour(np.flipud(possible[zz]), [0.5], colors="white", linewidths=0.5, extent=(0, zoom_km, 0, zoom_km))
    ax[4].set_title(f"Zoom {zoom_km} x {zoom_km} km (red box): radar groups underneath\n"
                    "red outline: burned, white: possibly burned", fontsize=10)
    ax[5].axis("off")
    rice_row = table[table["where"] == "rice"].iloc[0]
    other_row = table[table["where"] == "cropland, not rice"].iloc[0]
    ax[5].text(0, 0.95, "How to read this\n\n"
               "Each group is a set of cropland pixels whose radar curves look alike\n"
               "from June to October. Rice groups (greens) were dark in radar when\n"
               "flooded for transplanting and turned bright as the rice grew.\n\n"
               f"VIIRS fire points missed {rice_row.viirs_missed_share:.0%} of the burned area on rice\n"
               f"and {other_row.viirs_missed_share:.0%} on other cropland (step 2: 82% overall).", va="top", fontsize=11)
    for a in (ax[1], ax[4]):
        a.set_xlabel("km"); a.set_ylabel("km")
    fig.suptitle("Rice fields from Sentinel-1 radar, and the burn scars on them: rice belt north-west of Lahore", fontsize=13)
    fig.savefig(out, dpi=80)


if __name__ == "__main__":
    main()
