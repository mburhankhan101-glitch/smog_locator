"""
Step 5: Sheikhupura and Nankana Sahib districts, whole: what satellite fire points saw, what the burn
scars show, how much farmland burned, and how it lines up with smoke from Sentinel-5P (smoke.py).

With the environment variable SMOG_AREA=districts:

    python season.py run 2026; python season.py run 2025    # burn scars, both seasons (step 4's code)
    python districts.py fires                               # VIIRS fire points, both seasons (FIRMS_MAP_KEY)
    python smoke.py 2026-10-01 2026-11-01; python smoke.py 2025-10-01 2025-12-31
    python districts.py report                              # tables and figures in runs/districts/

Farmland is ESA WorldCover 2021 cropland. Districts are the union of their tehsils (geoBoundaries).

When fields burned. A scar is dated to a window: after the last clear Sentinel-2 look (about 11:00
local) up to the look that showed char. Fields burn in the afternoon, so each burned pixel is spread
evenly over the afternoons in its window: an estimate of hectares burned per day.

Fires against smoke. Burning, smoke and cold, calm nights all build up through October and November,
so any two of them look related. Besides the plain correlation, each daily series is also taken as its
difference from its own TREND_DAYS running median: on days with more fire than usual for the time of
the season, was there more smoke than usual?
"""
import os
import sys
from datetime import date as Date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage, stats

import aod
import burn_data
import burns
import fires
import season
import smoke

SEASONS = (2025, 2026)
TREND_DAYS = 15
COLOUR = {"Sheikhupura": "#2c3e50", "Nankana Sahib": "#16a085"}
SMOKE = {"aerosol_index": "UV aerosol index (Sentinel-5P)", "co": "carbon monoxide column (mmol/m², Sentinel-5P)"}


def season_end(year):
    return str(min(Date(year, 12, 31), Date.today()))


def district_grid():
    """District index of every pixel (-1: neither) and the district names."""
    from pyproj import Transformer
    from rasterio.features import rasterize
    from shapely.ops import transform
    tr, shape = burn_data.grid()
    to_utm = Transformer.from_crs("EPSG:4326", burn_data.CRS, always_xy=True).transform
    polygons = aod.districts()
    shapes = [(transform(to_utm, p), k) for k, p in enumerate(polygons.values())]
    return rasterize(shapes, out_shape=shape, transform=tr, fill=-1, dtype="int16"), list(polygons)


def one_season(year, district):
    """Burn scars, VIIRS fire points (persistent hot spots left out) and which scars a fire point saw."""
    s, rows, days = season.load(year)
    f = fires.on_grid(fires.viirs(f"{year}-10-01", season_end(year)))
    if len(f):
        f = f[~fires.persistent(f)].reset_index(drop=True)
    f["district"] = district[f.row, f.col]
    f["day"] = (f.time + fires.LOCAL).dt.strftime("%Y-%m-%d")
    b = s["burned_on"]
    hit = compared = np.zeros(b.shape, dtype=bool)
    if len(f):
        hit, _, lo, hi, _ = fires.compare(f, b, s["since"], days)
        first, last = fires.coverage((f.time.min(), f.time.max()), days)
        compared = (b >= 0) & (lo >= first) & (hi <= last)
        del lo, hi
    return {"year": year, "days": days, "table": pd.DataFrame(rows), "burned_on": b, "possible_on": s["possible_on"],
            "since": s["since"], "fires": f, "hit": hit & compared, "compared": compared}


def last_index(y, month_day):
    return max([k for k, d in enumerate(y["days"]) if d[5:] <= month_day], default=-1)


def district_table(ys, latest, district, names, crop):
    """One row per season, cut-off date and district (plus both together)."""
    rows = []
    areas = [(k, name) for k, name in enumerate(names)] + [(None, "both districts")]
    for y in ys:
        for up_to in sorted({latest, y["days"][-1][5:]} | ({"10-31"} if latest > "10-31" else set())):
            k = last_index(y, up_to)
            b = (y["burned_on"] >= 0) & (y["burned_on"] <= k)
            p = (y["possible_on"] >= 0) & (y["possible_on"] <= k)
            f = y["fires"][y["fires"].day.str[5:] <= up_to]
            for d, name in areas:
                m = district >= 0 if d is None else district == d
                compared = y["compared"] & b & m & crop                 # farmland only, as burned_cropland_ha
                rows.append({"season": y["year"], "up_to": up_to, "district": name,
                             "cropland_ha": (m & crop).sum() * burns.HA,
                             "viirs_fire_points": int((f.district >= 0).sum() if d is None else (f.district == d).sum()),
                             "burned_ha": (b & m).sum() * burns.HA,
                             "burned_cropland_ha": (b & m & crop).sum() * burns.HA,
                             "share_of_cropland_burned": (b & m & crop).sum() / max((m & crop).sum(), 1),
                             "possibly_burned_ha": (p & m).sum() * burns.HA,
                             "viirs_compared_ha": compared.sum() * burns.HA,
                             "viirs_caught_ha": (y["hit"] & compared).sum() * burns.HA,
                             "viirs_missed_share": 1 - (y["hit"] & compared).sum() / compared.sum()
                                                   if compared.any() else np.nan})
    return pd.DataFrame(rows)


def patches(ys, district, crop, bins=((1, 1), (2, 8), (9, 24), (25, 99), (100, None))):
    """Burned farmland by the size of its patch (touching burned pixels, diagonals included), and how much
    of each size VIIRS missed: is the result carried by tiny, uncertain patches?"""
    rows = []
    for y in ys:
        b = (y["burned_on"] >= 0) & (district >= 0) & crop
        label, _ = ndimage.label(b, structure=np.ones((3, 3)))
        size = np.bincount(label.ravel())[label]
        for lo, hi in bins:
            m = b & (size >= lo) & (size <= (hi or size.max()))
            compared = m & y["compared"]
            rows.append({"season": y["year"], "up_to": y["days"][-1][5:],
                         "patch_pixels": f"{lo}" if lo == hi else f"{lo}+" if hi is None else f"{lo}-{hi}",
                         "patch_ha": f"{lo * burns.HA:.2f}" + ("+" if hi is None else f"-{hi * burns.HA:.2f}"),
                         "burned_cropland_ha": m.sum() * burns.HA, "share_of_burned": m.sum() / b.sum(),
                         "viirs_compared_ha": compared.sum() * burns.HA,
                         "viirs_missed_share": 1 - (compared & y["hit"]).sum() / max(compared.sum(), 1)})
    return pd.DataFrame(rows)


def daily(y, district, names, crop, haze):
    """Per day and district: farmland burned (estimated, see above), VIIRS fire points, aerosol index, CO."""
    start = Date(y["year"], 10, 1)
    cal = pd.date_range(start, y["days"][-1]).strftime("%Y-%m-%d")
    since_day = np.array([Date.fromisoformat(d).toordinal() - start.toordinal() for d in y["days"]])
    out = []
    for d, name in enumerate(names):
        m = (district == d) & crop & (y["burned_on"] >= 0)
        pairs, counts = np.unique(np.stack([y["since"][m], y["burned_on"][m]]), axis=1, return_counts=True)
        burned = np.zeros(len(cal))
        for (s, k), c in zip(pairs.T, counts):
            a, z = since_day[s], since_day[k]               # afternoons of day a up to day z - 1
            burned[a:z] += c * burns.HA / (z - a)
        f = y["fires"][y["fires"].district == d].groupby("day").size()
        h = haze[haze.district == name].set_index("day")
        out.append(pd.DataFrame({"season": y["year"], "day": cal, "district": name, "burned_cropland_ha": burned,
                                 "fire_points": f.reindex(cal, fill_value=0).values,
                                 "aerosol_index": h.aerosol_index.reindex(cal).values, "co": h.co.reindex(cal).values}))
    return pd.concat(out, ignore_index=True)


def anomaly(x):
    return x - x.rolling(TREND_DAYS, center=True, min_periods=5).median()


def correlations(per_day):
    """Spearman correlation of the aerosol index and CO with fire points and with farmland burned, per
    season and district: plain, against the trend, and with the smoke one day later."""
    rows = []
    for (year, name), g in per_day.groupby(["season", "district"]):
        g = g.reset_index(drop=True)
        for smoke_series in ("aerosol_index", "co"):
            for what in ("fire_points", "burned_cropland_ha"):
                a = g[smoke_series]
                for how, x, s in (("plain", g[what], a), ("against the trend", anomaly(g[what]), anomaly(a)),
                                  ("smoke next day", g[what], a.shift(-1))):
                    ok = x.notna() & s.notna()
                    r, p = stats.spearmanr(x[ok], s[ok]) if ok.sum() > 5 else (np.nan, np.nan)
                    rows.append({"season": year, "district": name, "smoke": smoke_series, "series": what, "how": how,
                                 "days": int(ok.sum()), "spearman_r": r, "p_value": p})
    return pd.DataFrame(rows)


def weekly(per_day):
    w = per_day.copy()
    w["week"] = (pd.to_datetime(w.day) - pd.to_datetime(w.season.astype(str) + "-10-01")).dt.days // 7
    return (w.groupby(["season", "district", "week"])
             .agg(starts=("day", "first"), burned_cropland_ha=("burned_cropland_ha", "sum"),
                  fire_points=("fire_points", "sum"), aerosol_index=("aerosol_index", "mean"), co=("co", "mean"),
                  days=("day", "count"), smoke_days=("aerosol_index", "count"))
             .reset_index())


def report():
    import rice
    district, names = district_grid()
    crop = rice.worldcover() == rice.CROPLAND
    haze = smoke.load()
    ys = [one_season(year, district) for year in SEASONS if os.path.exists(season.STATE.format(year))]
    latest = min(y["days"][-1][5:] for y in ys)
    table = district_table(ys, latest, district, names, crop)
    per_day = pd.concat([daily(y, district, names, crop, haze) for y in ys], ignore_index=True)
    corr, per_week, sizes = correlations(per_day), weekly(per_day), patches(ys, district, crop)
    os.makedirs(burns.OUT, exist_ok=True)
    out = os.path.join(burns.OUT, "districts")
    for name, t in (("", table), ("_daily", per_day), ("_weekly", per_week), ("_correlation", corr),
                    ("_patches", sizes)):
        t.round(4).to_csv(out + name + ".csv", index=False)
    map_figure(ys, district, names, crop, out + "_map.png")
    bars_figure(table, latest, out + "_bars.png")
    series_figure(per_day, per_week, names, out + "_series.png")
    for col in SMOKE:
        scatter_figure(per_day, corr, out + f"_smoke_{col}.png", col)
    for name, t in (("districts", table), ("correlation", corr), ("patch sizes", sizes)):
        print(f"\n{name}\n" + t.round(3).to_string(index=False))


def map_figure(ys, district, names, crop, out, block=5):
    """Burned farmland (100 m blocks with any burned pixel) and VIIRS fire points, each season."""
    tr, shape = burn_data.grid()
    h, w = shape[0] // block * block, shape[1] // block * block
    pool = lambda a: a[:h, :w].reshape(h // block, block, w // block, block).any((1, 3))
    ext = (0, w * burn_data.RES / 1000, 0, h * burn_data.RES / 1000)
    fig, ax = plt.subplots(1, len(ys), figsize=(9 * len(ys), 9), layout="constrained", squeeze=False)
    for a, y in zip(ax[0], ys):
        a.imshow(np.where(pool(crop & (district >= 0)), 0.85, 1.0), cmap="gray", vmin=0, vmax=1, extent=ext)
        a.imshow(np.where(pool((y["burned_on"] >= 0) & crop & (district >= 0)), 1.0, np.nan), cmap="Reds", vmin=0,
                 vmax=1.1, extent=ext, interpolation="nearest")
        f = y["fires"][(y["fires"].district >= 0) & (y["fires"].day <= y["days"][-1])]     # same days as the scars
        a.scatter(f.col * burn_data.RES / 1000, ext[3] - f.row * burn_data.RES / 1000, s=5, c="orange",
                  edgecolors="none", alpha=0.9)
        a.contour(np.flipud(district[:h, :w][::block, ::block]), [-0.5, 0.5], colors="black", linewidths=1,
                  extent=ext)
        for k, name in enumerate(names):
            rr, cc = np.nonzero(district[::block, ::block] == k)
            a.text(np.median(cc) * block * burn_data.RES / 1000, ext[3] - rr.min() * block * burn_data.RES / 1000 - 4,
                   name, ha="center", fontsize=12, weight="bold")
        a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#cb181d", label="burned farmland (burn scars)"),
                          plt.Line2D([], [], color="orange", marker="o", ls="", label=f"VIIRS fire points ({len(f):,})"),
                          plt.Rectangle((0, 0), 1, 1, color="#d9d9d9", label="farmland (WorldCover)")],
                 loc="lower left", fontsize=10)
        a.set_title(f"{y['year']}: 1 October to {Date.fromisoformat(y['days'][-1]):%d %B}", fontsize=13)
        a.set_xlabel("km"); a.set_ylabel("km")
    fig.suptitle("Crop residue burning in Sheikhupura and Nankana Sahib: burn scars and fire points", fontsize=14)
    fig.savefig(out, dpi=80)


def bars_figure(table, latest, out):
    t = table[(table.up_to == latest) & (table.district != "both districts")]
    names, seasons = list(dict.fromkeys(t.district)), sorted(t.season.unique())
    fig, ax = plt.subplots(1, 2, figsize=(16, 6.5), layout="constrained")
    x, w = np.arange(len(names)), 0.38
    for i, year in enumerate(seasons):
        g = t[t.season == year].set_index("district").loc[names]
        ax[0].bar(x + (i - 0.5) * w, g.burned_cropland_ha, w, color="#cb181d", alpha=0.45 + 0.55 * i,
                  label=f"{year}: farmland burned (burn scars)")
        ax[0].bar(x + (i - 0.5) * w, g.viirs_caught_ha * g.burned_cropland_ha / g.viirs_compared_ha.clip(lower=1), w,
                  color="orange", alpha=0.45 + 0.55 * i, label=f"{year}: of which a VIIRS fire point saw")
        for xi, r in zip(x + (i - 0.5) * w, g.itertuples()):
            ax[0].text(xi, r.burned_cropland_ha, f"{r.burned_cropland_ha:,.0f} ha\n{r.share_of_cropland_burned:.1%} of farmland\n"
                       f"VIIRS missed {r.viirs_missed_share:.0%}", ha="center", va="bottom", fontsize=9)
        ax[1].bar(x + (i - 0.5) * w, g.viirs_fire_points, w, color="orange", alpha=0.45 + 0.55 * i, label=str(year))
    ax[0].set_xticks(x, names, fontsize=11); ax[1].set_xticks(x, names, fontsize=11)
    ax[0].set_ylim(0, t.burned_cropland_ha.max() * 1.3)
    ax[0].set_ylabel("hectares"); ax[1].set_ylabel("VIIRS fire points")
    ax[0].legend(fontsize=9, loc="upper right"); ax[1].legend(fontsize=9)
    day = f"{int(latest[3:])} {pd.Timestamp(2000, int(latest[:2]), 1):%B}"
    ax[0].set_title(f"Farmland burned, 1 October to {day}, and how much fire points caught", fontsize=12)
    ax[1].set_title(f"VIIRS fire points, 1 October to {day}", fontsize=12)
    fig.savefig(out, dpi=80)


def series_figure(per_day, per_week, names, out):
    seasons = sorted(per_day.season.unique())
    fig, ax = plt.subplots(4, len(seasons), figsize=(10 * len(seasons), 15), layout="constrained", squeeze=False,
                           sharey="row")
    for j, year in enumerate(seasons):
        d = per_day[per_day.season == year]
        both = d.groupby("day")[["burned_cropland_ha", "fire_points"]].sum()
        x = pd.to_datetime(both.index)
        ax[0, j].bar(x, both.burned_cropland_ha, color="#cb181d", width=0.9, label="farmland burned (estimated, ha)")
        tw = ax[0, j].twinx()
        tw.plot(x, both.fire_points, color="orange", marker="o", ms=3, lw=1.5, label="VIIRS fire points")
        tw.set_ylabel("VIIRS fire points per day", color="darkorange")
        ax[0, j].set_ylabel("farmland burned per day (ha)")
        ax[0, j].set_title(f"{year}: both districts, per day", fontsize=12)
        h1, l1 = ax[0, j].get_legend_handles_labels(); h2, l2 = tw.get_legend_handles_labels()
        ax[0, j].legend(h1 + h2, l1 + l2, fontsize=9, loc="upper left")
        for i, (col, label) in enumerate(SMOKE.items(), start=1):
            for name in names:
                g = d[d.district == name]
                gx = pd.to_datetime(g.day)
                ax[i, j].scatter(gx, g[col], s=14, color=COLOUR[name])
                ax[i, j].plot(gx, g[col].interpolate(limit=3).rolling(5, center=True, min_periods=2).mean(),
                              color=COLOUR[name], lw=1.5, label=name)
            ax[i, j].set_ylabel(label)
            ax[i, j].legend(fontsize=9, loc="upper left")
            ax[i, j].set_title(f"{label.split(' (')[0].capitalize()}, district average per day at about 13:30 "
                               "(dots; line: 5-day mean; gaps: clouds)", fontsize=11)
        wk = per_week[per_week.season == year]
        wb = wk.groupby("week").agg(starts=("starts", "first"), burned=("burned_cropland_ha", "sum"),
                                    days=("days", "first"), aerosol_index=("aerosol_index", "mean"), co=("co", "mean"))
        wx = pd.to_datetime(wb.starts) + pd.Timedelta(days=3)
        full = wb.days == 7
        ax[3, j].bar(wx[full], wb.burned[full], width=6, color="#cb181d", alpha=0.8, label="farmland burned (ha)")
        if (~full).any():
            ax[3, j].bar(wx[~full], wb.burned[~full], width=6, color="#cb181d", alpha=0.3, hatch="//",
                         label="week not over yet")
        for xi, r in zip(wx[~full], wb[~full].itertuples()):
            ax[3, j].text(xi, r.burned, f"{r.days} of 7 days", ha="center", va="bottom", fontsize=9)
        tw = ax[3, j].twinx()
        tw.plot(wx, wb.aerosol_index, color="black", marker="s", lw=2, label="aerosol index, weekly mean")
        tw.set_ylabel("aerosol index")
        tc = ax[3, j].twinx()
        tc.spines.right.set_position(("axes", 1.07))
        tc.plot(wx, wb.co, color="grey", marker="o", ls="--", lw=1.5, label="CO, weekly mean")
        tc.set_ylabel("CO (mmol/m²)", color="dimgrey")
        ax[3, j].set_ylabel("farmland burned per week (ha)")
        ax[3, j].set_title("Per week (both districts)", fontsize=12)
        hs = [h for a in (ax[3, j], tw, tc) for h in a.get_legend_handles_labels()[0]]
        ax[3, j].legend(hs, [h.get_label() for h in hs], fontsize=9, loc="upper left")
        for a in ax[:, j]:
            a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
            a.set_xlim(pd.Timestamp(year, 9, 30), pd.Timestamp(year, 12, 31))
    fig.suptitle("Burning and smoke through the season, Sheikhupura and Nankana Sahib", fontsize=14)
    fig.savefig(out, dpi=80)


def scatter_figure(per_day, corr, out, col="aerosol_index"):
    seasons = sorted(per_day.season.unique())
    fig, ax = plt.subplots(2, len(seasons), figsize=(8 * len(seasons), 12), layout="constrained", squeeze=False)
    for j, year in enumerate(seasons):
        for name, g in per_day[per_day.season == year].groupby("district"):
            g = g.reset_index(drop=True)
            for i, (x, a, how) in enumerate(((g.fire_points, g[col], "plain"),
                                             (anomaly(g.fire_points), anomaly(g[col]), "against the trend"))):
                c = corr[(corr.season == year) & (corr.district == name) & (corr.smoke == col)
                         & (corr.series == "fire_points") & (corr.how == how)].iloc[0]
                ax[i, j].scatter(x, a, s=18, color=COLOUR[name], alpha=0.8,
                                 label=f"{name}: r = {c.spearman_r:.2f} ({c.days} days, p = {c.p_value:.2f})")
        short = SMOKE[col].split(" (")[0]
        ax[0, j].set_xlabel("VIIRS fire points that day"); ax[0, j].set_ylabel(f"{short} that day")
        ax[1, j].set_xlabel(f"fire points minus their {TREND_DAYS}-day running median")
        ax[1, j].set_ylabel(f"{short} minus its {TREND_DAYS}-day running median")
        ax[0, j].set_title(f"{year}: fire points against smoke, day by day", fontsize=12)
        ax[1, j].set_title(f"{year}: the same with the seasonal build-up taken out\n"
                           "(more fire than usual for the time of season: more smoke than usual?)", fontsize=12)
        for a in ax[:, j]:
            a.legend(fontsize=9, loc="upper left")
            a.axvline(0, color="grey", lw=0.5)
        ax[1, j].axhline(0, color="grey", lw=0.5)
    fig.suptitle(f"Do more fires go with more smoke? Spearman correlation per district ({SMOKE[col]})", fontsize=14)
    fig.savefig(out, dpi=80)


def fetch_fires():
    for year in SEASONS:
        f = fires.viirs(f"{year}-10-01", season_end(year))
        print(f"{year}: {len(f):,} VIIRS fire points over the area, {f.time.min()} to {f.time.max()}", flush=True)


if __name__ == "__main__":
    if burn_data.AREA != "districts":
        sys.exit("set SMOG_AREA=districts first")
    {"fires": fetch_fires, "report": report}[sys.argv[1]]()
