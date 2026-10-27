"""
Satellite fire points against the burn scars of step 1: how much burning do fire counts miss?

Fire points are hot spots seen in the middle-infrared (about 4 um, where a flame at 800-1200 K is
far brighter than the ground at 300 K) at the moment a satellite passes:

    VIIRS        375 m pixels, on three satellites (Suomi NPP, NOAA-20, NOAA-21), passing at about
                 13:30 and 01:30 local time; from NASA FIRMS. The public files only hold the last 7
                 days, so every download is kept in cache/firms/ and they are merged. With a FIRMS
                 MAP_KEY in the environment variable FIRMS_MAP_KEY, older days are fetched too.
    Sentinel-3   SLSTR fire product, 1 km pixels, passing at about 10:00 and 22:00 local time; from
                 Microsoft Planetary Computer. Only "vegetation fire" points are kept.

Persistent hot spots, mostly brick kilns, fire on many days in the same place and leave no burn
scar: a 400 m cell with fire points on PERSISTENT_DAYS or more different days is left out.

A burned pixel (step 1) burned between its last clear Sentinel-2 look and the day it showed char,
both at about S2_UTC. It counts as "seen by a fire point" when a fire point fell within the fire
pixel's size (RADIUS) of it in that window. Only burned pixels whose whole window lies inside the
days a sensor covers are compared. The other way round, a fire point "has a burn scar" when a
pixel within RADIUS was flagged in the window around the fire's time.

    python fires.py 2026-10-01 2026-10-27        # needs runs/burns_<start>_<end>.npz from burns.py
"""
import glob
import io
import os
import sys
from datetime import date as Date, datetime, timedelta, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from pyproj import Transformer

import burn_data
import burns

FIRMS = "https://firms.modaps.eosdis.nasa.gov"
VIIRS_7D = {"N": "suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_South_Asia_7d.csv",
            "N20": "noaa-20-viirs-c2/csv/J1_VIIRS_C2_South_Asia_7d.csv",
            "N21": "noaa-21-viirs-c2/csv/J2_VIIRS_C2_South_Asia_7d.csv"}
VIIRS_API = ("VIIRS_SNPP_NRT", "VIIRS_NOAA20_NRT", "VIIRS_NOAA21_NRT")
RADIUS = {"VIIRS": 375, "S3": 1000}
S2_UTC = 6                      # Sentinel-2 passes the area at about 05:55 UTC (10:55 local)
LOCAL = timedelta(hours=5)      # Pakistan time
PERSISTENT_DAYS, CELL = 3, 400


def viirs():
    """VIIRS fire points over the area: every saved 7-day file, plus the archive with a MAP_KEY."""
    folder = os.path.join(burn_data.CACHE, "firms")
    os.makedirs(folder, exist_ok=True)
    today = str(Date.today())
    for sat, path in VIIRS_7D.items():
        out = os.path.join(folder, f"{sat}_7d_{today}.csv")
        if not os.path.exists(out):
            r = requests.get(f"{FIRMS}/data/active_fire/{path}", timeout=300)
            r.raise_for_status()
            open(out, "wb").write(r.content)
    key = os.environ.get("FIRMS_MAP_KEY")
    if key:
        area = ",".join(str(v) for v in burn_data.LONLAT)
        d = Date(2026, 10, 1)
        while d < Date.today():
            for src in VIIRS_API:
                out = os.path.join(folder, f"{src}_{d}.csv")
                if not os.path.exists(out):
                    r = requests.get(f"{FIRMS}/api/area/csv/{key}/{src}/{area}/5/{d}", timeout=300)
                    r.raise_for_status()
                    open(out, "wb").write(r.content)
            d += timedelta(days=5)
    w, s, e, n = burn_data.LONLAT
    f = pd.concat([pd.read_csv(p, dtype={"acq_time": str}) for p in glob.glob(os.path.join(folder, "*.csv"))])
    f = f[(f.longitude > w) & (f.longitude < e) & (f.latitude > s) & (f.latitude < n)]
    f = f[f.confidence.astype(str).str.lower().str[0] != "l"]          # drop low confidence
    f = f.drop_duplicates(["latitude", "longitude", "acq_date", "acq_time", "satellite"])
    t = pd.to_datetime(f.acq_date + " " + f.acq_time.str.zfill(4), format="%Y-%m-%d %H%M", utc=True)
    return pd.DataFrame({"sensor": "VIIRS", "satellite": f.satellite.astype(str).values, "time": t.values,
                         "lat": f.latitude.values, "lon": f.longitude.values, "frp": f.frp.values})


def s3_items(start, end, chunk=3, tries=4):
    """Sentinel-3 fire product passes over the area, searched a few days at a time (long searches time out)."""
    import pystac_client
    d, last = Date.fromisoformat(start), Date.fromisoformat(end)
    while d <= last:
        e = min(d + timedelta(days=chunk - 1), last)
        for k in range(tries):
            try:
                items = list(burn_data.catalog().search(collections=["sentinel-3-slstr-frp-l2-netcdf"],
                                                        bbox=burn_data.LONLAT, datetime=f"{d}/{e}T23:59:59Z",
                                                        method="GET", limit=100).items())
                break
            except pystac_client.exceptions.APIError:
                if k == tries - 1:
                    raise
        yield from items
        d = e + timedelta(days=1)


def sentinel3(start, end):
    """Sentinel-3 vegetation fire points over the area, one small CSV per pass in cache/s3frp/,
    and the times of the first and last pass."""
    import xarray as xr
    folder = os.path.join(burn_data.CACHE, "s3frp")
    os.makedirs(folder, exist_ok=True)
    w, s, e, n = burn_data.LONLAT
    parts, passes = [], []
    for item in s3_items(start, end):
        passes.append(pd.Timestamp(item.datetime))
        out = os.path.join(folder, item.id + ".csv")
        if not os.path.exists(out):
            r = requests.get(item.assets["frp-in"].href, timeout=300)
            r.raise_for_status()
            ds = xr.open_dataset(io.BytesIO(r.content), engine="h5netcdf")
            keep = ((ds.longitude > w) & (ds.longitude < e) & (ds.latitude > s) & (ds.latitude < n)
                    & (ds.classification == 1)).values
            pd.DataFrame({"time": ds.time.values[keep], "lat": ds.latitude.values[keep],
                          "lon": ds.longitude.values[keep], "frp": ds.FRP_MWIR.values[keep]}).to_csv(out, index=False)
        p = pd.read_csv(out, parse_dates=["time"])
        p["satellite"] = item.id[:3]
        parts.append(p)
    f = pd.concat(parts)
    f["time"] = pd.to_datetime(f.time, utc=True)
    f.insert(0, "sensor", "S3")
    return f.reset_index(drop=True), (min(passes), max(passes))


def on_grid(f):
    """Row and column of each fire point on the 20 m grid."""
    tr, shape = burn_data.grid()
    x, y = Transformer.from_crs("EPSG:4326", burn_data.CRS, always_xy=True).transform(f.lon.values, f.lat.values)
    f = f.copy()
    f["time"] = pd.to_datetime(f.time, utc=True)
    f["x"], f["y"] = x, y
    f["row"] = ((tr.f - y) / burn_data.RES).astype(int)
    f["col"] = ((x - tr.c) / burn_data.RES).astype(int)
    return f[(f.row >= 0) & (f.row < shape[0]) & (f.col >= 0) & (f.col < shape[1])].reset_index(drop=True)


def persistent(f):
    """True for fire points in a CELL x CELL m cell that had fire points on PERSISTENT_DAYS or more days."""
    cell = (f.x // CELL).astype(int).astype(str) + "_" + (f.y // CELL).astype(int).astype(str)
    days = f.groupby(cell)["time"].transform(lambda t: t.dt.date.nunique())
    return (days >= PERSISTENT_DAYS).values


def disc(row, col, radius, shape):
    """Slices and mask of the pixels within radius (m) of a fire point."""
    r = int(np.ceil(radius / burn_data.RES))
    r0, r1, c0, c1 = max(row - r, 0), min(row + r + 1, shape[0]), max(col - r, 0), min(col + r + 1, shape[1])
    yy, xx = np.mgrid[r0:r1, c0:c1]
    return (slice(r0, r1), slice(c0, c1)), (yy - row) ** 2 + (xx - col) ** 2 <= r * r


def compare(f, flagged_on, since, days):
    """Which flagged pixels had a fire point in their burn window, and which fire points had a scar.
    Times are hours since the start; a pixel's window is (last clear look, day flagged]."""
    t0 = datetime.fromisoformat(days[0]).replace(tzinfo=timezone.utc)
    s2 = np.array([(datetime.fromisoformat(d).replace(tzinfo=timezone.utc) - t0).total_seconds() / 3600 + S2_UTC
                   for d in days])
    hit = np.zeros(flagged_on.shape, dtype=bool)
    scar = np.zeros(len(f), dtype=bool)
    lo = np.where(flagged_on >= 0, s2[np.maximum(since, 0)], np.inf)
    hi = np.where(flagged_on >= 0, s2[np.maximum(flagged_on, 0)], -np.inf)
    hours = ((f.time - pd.Timestamp(t0)).dt.total_seconds() / 3600).values
    for k, (row, col, sensor, h) in enumerate(zip(f.row, f.col, f.sensor, hours)):
        z, m = disc(row, col, RADIUS[sensor], flagged_on.shape)
        inside = m & (lo[z] < h) & (hi[z] >= h)
        hit[z] |= inside
        scar[k] = inside.any()
    return hit, scar, lo, hi, hours


def coverage(span, days):
    """First and last time a sensor covers, in hours since the start."""
    t0 = pd.Timestamp(days[0], tz="UTC")
    return tuple((pd.Timestamp(t).tz_convert("UTC") - t0).total_seconds() / 3600 for t in span)


def main():
    start, end = sys.argv[1], sys.argv[2]
    z = np.load(os.path.join(burns.OUT, f"burns_{start}_{end}.npz"))
    burned_on, possible_on, since, days = z["burned_on"], z["possible_on"], z["since"], list(z["days"])
    v = viirs()
    s3, s3_span = sentinel3(start, end)
    v_time = pd.to_datetime(v.time, utc=True)
    span = {"VIIRS": (v_time.min(), v_time.max()), "S3": s3_span}     # VIIRS: its first and last point
    f = on_grid(pd.concat([v, s3], ignore_index=True))
    f["persistent"] = persistent(f)
    local = f.time + LOCAL
    f["local_hour"] = (local.dt.hour + local.dt.minute / 60).values
    s2_last = (Date.fromisoformat(days[-1]) - Date.fromisoformat(days[0])).days * 24 + S2_UTC
    rows, hits, compared = [], {}, {}
    for sensor in ("VIIRS", "S3"):
        g = f[(f.sensor == sensor) & ~f.persistent].reset_index(drop=True)
        first, last = coverage(span[sensor], days)
        scar = np.zeros(len(g), dtype=bool)
        for level, flagged in (("burned", burned_on), ("possibly burned", possible_on)):
            hit, scar_k, lo, hi, hours = compare(g, flagged, since, days)
            inside = (flagged >= 0) & (lo >= first) & (hi <= last)
            hits[sensor, level], compared[sensor, level] = hit & inside, inside
            scar |= scar_k
            rows.append({"sensor": sensor, "level": level,
                         "covered": f"{(pd.Timestamp(days[0]) + pd.Timedelta(hours=first)):%d %b} to "
                                    f"{(pd.Timestamp(days[0]) + pd.Timedelta(hours=last)):%d %b}",
                         "flagged_ha_compared": inside.sum() * burns.HA,
                         "seen_by_fire_point_ha": (hit & inside).sum() * burns.HA,
                         "missed_share": 1 - (hit & inside).sum() / max(inside.sum(), 1)})
        in_s2 = (hours > S2_UTC) & (hours <= s2_last)       # fire points between the first and last S2 look
        for r in rows[-2:]:
            r["fire_points"], r["fire_points_with_scar"] = int(in_s2.sum()), int((scar & in_s2).sum())
    table = pd.DataFrame(rows)
    os.makedirs(burns.OUT, exist_ok=True)
    stem = os.path.join(burns.OUT, f"fires_{start}_{end}")
    f.drop(columns=["x", "y"]).to_csv(stem + "_points.csv", index=False)
    table.round(3).to_csv(stem + ".csv", index=False)
    figure(days, burned_on, possible_on, f, table, hits, compared, stem + ".png")
    print(f"fire points: {len(f)} ({(f.sensor == 'VIIRS').sum()} VIIRS, {(f.sensor == 'S3').sum()} Sentinel-3), "
          f"{f.persistent.sum()} at persistent hot spots")
    print(table.round(3).to_string(index=False))


def figure(days, burned_on, possible_on, f, table, hits, compared, out, zoom_km=4):
    tr, shape = burn_data.grid()
    w_km, h_km = shape[1] * burn_data.RES / 1000, shape[0] * burn_data.RES / 1000
    km = lambda g: (g.col * burn_data.RES / 1000, h_km - g.row * burn_data.RES / 1000)
    last = burn_data.day(days[-1])[0]
    fig = plt.figure(figsize=(20, 13), layout="constrained")
    g = fig.add_gridspec(2, 3)
    ax = [fig.add_subplot(g[0, 0]), fig.add_subplot(g[0, 1]), fig.add_subplot(g[0, 2]),
          fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1:])]
    base = burns.false_colour(last)[::3, ::3].mean(-1)
    ax[0].imshow(base, cmap="gray", vmin=0, vmax=1, extent=(0, w_km, 0, h_km))
    ax[0].imshow(np.where(possible_on >= 0, 1.0, np.nan), cmap="Greys", vmin=0, vmax=1.5, alpha=0.6,
                 extent=(0, w_km, 0, h_km), interpolation="nearest")
    ax[0].imshow(np.where(burned_on >= 0, 1.0, np.nan), cmap="Reds", vmin=0, vmax=1.2,
                 extent=(0, w_km, 0, h_km), interpolation="nearest")
    for sensor, colour, marker in (("VIIRS", "orange", "o"), ("S3", "deepskyblue", "s")):
        p = f[(f.sensor == sensor) & ~f.persistent]
        ax[0].scatter(*km(p), s=8, c=colour, marker=marker, label=f"{sensor} fire points ({len(p)})")
    p = f[f.persistent]
    ax[0].scatter(*km(p), s=30, c="black", marker="x", label=f"persistent hot spots, left out ({len(p)})")
    ax[0].legend(fontsize=9, loc="upper right")
    ax[0].set_title(f"Burned fields (red), possibly burned (grey) and\nsatellite fire points, {days[0]} to {days[-1]}",
                    fontsize=10)
    # share of compared burned area that a fire point saw
    labels, seen, missed = [], [], []
    for r in table.itertuples():
        labels.append(f"{r.sensor}\n{r.level}\n{r.covered}")
        seen.append(r.seen_by_fire_point_ha)
        missed.append(r.flagged_ha_compared - r.seen_by_fire_point_ha)
    x = np.arange(len(labels))
    ax[1].bar(x, seen, color="#e67e22", label="a fire point fell on it")
    ax[1].bar(x, missed, bottom=seen, color="#7f8c8d", label="no fire point: missed")
    for k, r in enumerate(table.itertuples()):
        ax[1].text(k, r.flagged_ha_compared, f"{r.missed_share:.0%}\nmissed", ha="center", va="bottom", fontsize=9)
    ax[1].set_xticks(x, labels, fontsize=8)
    ax[1].set_ylabel("flagged area compared (ha)")
    ax[1].legend(fontsize=8, loc="upper left")
    ax[1].set_title("How much of the burned area did fire points catch?\n(only fields burned while the sensor's data covers)",
                    fontsize=10)
    # local time of the fire points against the Sentinel-2 look
    bins = np.arange(0, 25)
    for sensor, colour in (("VIIRS", "orange"), ("S3", "deepskyblue")):
        ax[2].hist(f[f.sensor == sensor].local_hour, bins=bins, color=colour, alpha=0.8, label=sensor)
    ax[2].axvline(S2_UTC + 5, color="red", lw=1, ls="--", label="Sentinel-2 look (burn scars)")
    ax[2].axvspan(15, 19, color="grey", alpha=0.2, label="late afternoon: no fire satellite passes")
    ax[2].set_xticks(range(0, 25, 3))
    ax[2].set_xlabel("local time (hour)"); ax[2].set_ylabel("fire points")
    ax[2].legend(fontsize=8)
    ax[2].set_title("Fire points only exist at the moments a satellite passes", fontsize=10)
    # zoom: the window with the most burned area VIIRS could be compared with, fire point footprints drawn
    flagged, caught = compared["VIIRS", "burned"], hits["VIIRS", "burned"]
    n = int(zoom_km * 1000 / burn_data.RES)
    dens = np.cumsum(np.cumsum(np.pad(flagged.astype("int32"), ((1, 0), (1, 0))), 0), 1)
    box = dens[n:, n:] - dens[:-n, n:] - dens[n:, :-n] + dens[:-n, :-n]
    r0, c0 = np.unravel_index(np.argmax(box), box.shape)
    zz = (slice(r0, r0 + n), slice(c0, c0 + n))
    x0, y0 = c0 * burn_data.RES / 1000, h_km - (r0 + n) * burn_data.RES / 1000
    ax[0].add_patch(plt.Rectangle((x0, y0), zoom_km, zoom_km, fill=False, edgecolor="cyan", lw=1.5))
    ax[3].imshow(burns.false_colour(last)[zz], extent=(0, zoom_km, 0, zoom_km), interpolation="nearest")
    ax[3].contour(np.flipud(flagged[zz]), [0.5], colors="white", linewidths=0.6, extent=(0, zoom_km, 0, zoom_km))
    ax[3].contour(np.flipud(caught[zz]), [0.5], colors="orange", linewidths=1.2, extent=(0, zoom_km, 0, zoom_km))
    for sensor, colour in (("VIIRS", "orange"), ("S3", "deepskyblue")):
        p = f[(f.sensor == sensor) & ~f.persistent & (f.row >= r0) & (f.row < r0 + n) & (f.col >= c0) & (f.col < c0 + n)]
        for row, col in zip(p.row, p.col):
            ax[3].add_patch(plt.Circle(((col - c0) * burn_data.RES / 1000, zoom_km - (row - r0) * burn_data.RES / 1000),
                                       RADIUS[sensor] / 1000, fill=False, edgecolor=colour, lw=1))
    ax[3].set_xlim(0, zoom_km); ax[3].set_ylim(0, zoom_km)
    ax[3].set_title(f"Zoom {zoom_km} x {zoom_km} km (cyan box), {days[-1]}: fields burned 21-26 Oct\n"
                    "white: burned, orange: burned and caught by a VIIRS fire point (circle)", fontsize=10)
    ax[4].axis("off")
    ax[4].text(0, 0.95, "How to read this\n\n"
               "A fire point is a hot spot a satellite saw at the moment it passed:\n"
               "VIIRS (375 m) at about 13:30 and 01:30, Sentinel-3 (1 km) at about\n"
               "10:00 and 22:00. A field that burned at 16:00 and was out by 18:00\n"
               "gives no fire point at all.\n\n"
               "A burn scar stays black for days, whenever the field burned.\n"
               "Sentinel-2 sees it at its next clear look.\n\n"
               "Persistent hot spots (fire on 3 or more days in the same 400 m\n"
               "cell) are mostly brick kilns and are left out.", va="top", fontsize=11)
    for a in (ax[0], ax[3]):
        a.set_xlabel("km"); a.set_ylabel("km")
    fig.suptitle("What satellite fire counts miss: fire points against burn scars, rice belt north-west of Lahore",
                 fontsize=13)
    fig.savefig(out, dpi=80)


if __name__ == "__main__":
    main()
