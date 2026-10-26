"""
Burned fields from Sentinel-2: every pixel is compared with its own last clear look.

    NBR    (NIR - SWIR2) / (NIR + SWIR2), with the narrow NIR band B8A and B12; drops when a field burns
    MIRBI  10 x SWIR2 - 9.8 x SWIR1 + 2 (reflectance 0-1), an index built for crop and grass fires;
           rises when a field burns

Going through the days in order, each valid pixel keeps its last clear NBR, MIRBI, water flag and
day. Against its last clear look no more than MAX_GAP days earlier, a pixel is flagged on day d as

    burned           dNBR = NBR before - NBR now >= DNBR_MIN (a big drop), NBR now <= NBR_MAX and
                     NIR now <= NIR_MAX: very dark in NIR, char rather than harvested straw
    possibly burned  the same drop and NBR, NIR between NIR_MAX and NIR_POSSIBLE (char mixed with
                     bright ash, or something else), and MIRBI up by at least DMIRBI_MIN

Why these rules (checked on 18 to 21 October 2026): about 68,000 ha had a big NBR drop in those three
days, mostly harvest and haze; harvested straw stays fairly bright in NIR (median 0.27), char does
not. Images of the flagged fields show straw or bare fields turning dark maroon in false colour.
MIRBI rose by a median 0.67 on "burned" pixels and 0.20 on ordinary harvest. Water (MNDWI > 0 before
or after, plus 2 pixels around it) is left out: shifting river and canal edges looked like burns.

The burn happened between the last clear look and day d. Each pixel counts once (its first
detection). Area in hectares: 20 m pixel = 0.04 ha.

    python burns.py 2026-10-01 2026-10-27          # writes runs/burns_<start>_<end>.npz, .csv, .png
"""
import os
import sys
from datetime import date as Date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import burn_data

DNBR_MIN, NBR_MAX, NIR_MAX, NIR_POSSIBLE, DMIRBI_MIN, MAX_GAP = 0.20, 0.0, 0.15, 0.25, 0.4, 16
HA = (burn_data.RES ** 2) / 10_000
OUT = "runs"


def indices(refl):
    """NBR, NIR, MIRBI and a water flag (MNDWI > 0) for one day."""
    r = refl.astype("float32") / 10000
    green, nir, swir1, swir2 = r[1], r[3], r[4], r[5]
    with np.errstate(invalid="ignore", divide="ignore"):
        nbr = (nir - swir2) / (nir + swir2)
        water = (green - swir1) / (green + swir1) > 0
    return nbr, nir, 10 * swir2 - 9.8 * swir1 + 2, water


def detect(days, dnbr_min=DNBR_MIN, nbr_max=NBR_MAX, nir_max=NIR_MAX, max_gap=MAX_GAP):
    """burned_on / possible_on: index of the day each pixel was first flagged (-1: never);
    since: the last clear look before it; and a table per day."""
    from scipy import ndimage
    tr, shape = burn_data.grid()
    last_nbr = np.full(shape, np.nan, dtype="float32")
    last_mirbi = np.full(shape, np.nan, dtype="float32")
    last_water = np.zeros(shape, dtype=bool)
    possible_on = np.full(shape, -1, dtype="int16")
    last_day = np.full(shape, -1, dtype="int32")
    burned_on = np.full(shape, -1, dtype="int16")      # index of the day it was first seen burned
    since = np.full(shape, -1, dtype="int16")          # index of the last clear look before that
    rows = []
    day_index = {d: k for k, d in enumerate(days)}
    ordinal = {d: Date.fromisoformat(d).toordinal() for d in days}
    for d in days:
        refl, valid = burn_data.day(d)
        nbr, nir, mirbi, water = indices(refl)
        k = day_index[d]
        gap = ordinal[d] - np.where(last_day >= 0, np.array([ordinal[x] for x in days])[np.maximum(last_day, 0)], 0)
        seen = valid & (last_day >= 0) & (gap <= max_gap)
        wet = ndimage.binary_dilation((water & valid) | last_water, iterations=2)
        drop = seen & ~wet & (last_nbr - nbr >= dnbr_min) & (nbr <= nbr_max) & (burned_on < 0) & (possible_on < 0)
        new = drop & (nir <= nir_max)
        maybe = drop & (nir > nir_max) & (nir <= NIR_POSSIBLE) & (mirbi - last_mirbi >= DMIRBI_MIN)
        burned_on[new] = k
        possible_on[maybe] = k
        since[new | maybe] = last_day[new | maybe]
        rows.append({"day": d, "valid_share": float(valid.mean()), "compared_share": float(seen.mean()),
                     "newly_burned_ha": float(new.sum() * HA), "newly_possibly_burned_ha": float(maybe.sum() * HA)})
        last_nbr[valid] = nbr[valid]
        last_mirbi[valid] = mirbi[valid]
        last_water[valid] = water[valid]
        last_day[valid] = k
    return burned_on, possible_on, since, pd.DataFrame(rows)


def false_colour(refl):
    """SWIR-2, NIR, red as red, green, blue: char black, vegetation green, bare soil pink-brown."""
    r = refl.astype("float32") / 10000
    return np.stack([r[5] / 0.35, r[3] / 0.45, r[2] / 0.2], -1).clip(0, 1)


def figure(days, burned_on, possible_on, table, out, zoom_km=5):
    tr, shape = burn_data.grid()
    first, last = burn_data.day(days[0])[0], burn_data.day(days[-1])[0]
    w_km, h_km = shape[1] * burn_data.RES / 1000, shape[0] * burn_data.RES / 1000
    fig = plt.figure(figsize=(20, 13), layout="constrained")
    g = fig.add_gridspec(2, 3)
    ax = [fig.add_subplot(g[0, 0]), fig.add_subplot(g[0, 1]), fig.add_subplot(g[0, 2]),
          fig.add_subplot(g[1, 0]), fig.add_subplot(g[1, 1]), fig.add_subplot(g[1, 2])]
    ax[0].imshow(false_colour(last)[::3, ::3], extent=(0, w_km, 0, h_km))
    ax[0].set_title(f"False colour (SWIR-2, NIR, red), {days[-1]}\nchar is black or dark maroon, crops green, bare soil pink",
                    fontsize=10)
    cmap = plt.get_cmap("plasma", len(days))
    shown = np.where(burned_on >= 0, burned_on, np.nan).astype("float32")
    ax[1].imshow(false_colour(last)[::3, ::3].mean(-1), cmap="gray", vmin=0, vmax=1, extent=(0, w_km, 0, h_km))
    ax[1].imshow(np.where(possible_on >= 0, 1.0, np.nan), cmap="Greys", vmin=0, vmax=1.5, alpha=0.6,
                 extent=(0, w_km, 0, h_km), interpolation="nearest")
    im = ax[1].imshow(shown, cmap=cmap, vmin=-0.5, vmax=len(days) - 0.5, extent=(0, w_km, 0, h_km), interpolation="nearest")
    cb = fig.colorbar(im, ax=ax[1], shrink=0.75, ticks=range(len(days)))
    cb.ax.set_yticklabels(days, fontsize=7)
    ax[1].set_title(f"Burned fields, coloured by the day first seen (grey: possibly burned)\n"
                    f"{(burned_on >= 0).sum() * HA:,.0f} ha burned, {(possible_on >= 0).sum() * HA:,.0f} ha possibly",
                    fontsize=10)
    ax[2].bar(range(len(table)), table.newly_burned_ha, color="#c0392b", label="burned")
    ax[2].bar(range(len(table)), table.newly_possibly_burned_ha, bottom=table.newly_burned_ha, color="#bbbbbb",
              label="possibly burned")
    ax[2].legend(fontsize=8, loc="upper left")
    ax[2].set_xticks(range(len(table)), table.day, rotation=60, fontsize=8)
    ax[2].set_ylabel("newly flagged (ha)")
    ax2 = ax[2].twinx()
    ax2.plot(range(len(table)), table.compared_share * 100, color="grey", marker="o", lw=1)
    ax2.set_ylabel("share of the area compared (%), grey", color="grey")
    ax[2].set_title("Newly flagged area per Sentinel-2 day\n(low grey = cloudy day: little could be checked)", fontsize=10)
    # the 5 x 5 km window with the most burned area: first day, last day, flagged fields outlined
    n = int(zoom_km * 1000 / burn_data.RES)
    dens = np.cumsum(np.cumsum(np.pad((burned_on >= 0).astype("int32"), ((1, 0), (1, 0))), 0), 1)
    box = dens[n:, n:] - dens[:-n, n:] - dens[n:, :-n] + dens[:-n, :-n]
    r0, c0 = np.unravel_index(np.argmax(box), box.shape)
    z = (slice(r0, r0 + n), slice(c0, c0 + n))
    x0, y0 = c0 * burn_data.RES / 1000, h_km - (r0 + n) * burn_data.RES / 1000
    ax[1].add_patch(plt.Rectangle((x0, y0), zoom_km, zoom_km, fill=False, edgecolor="cyan", lw=1.5))
    for a, img, title in ((ax[3], first, days[0]), (ax[4], last, days[-1])):
        a.imshow(false_colour(img)[z], extent=(0, zoom_km, 0, zoom_km), interpolation="nearest")
        a.contour(np.flipud(burned_on[z] >= 0), [0.5], colors="cyan", linewidths=0.8, extent=(0, zoom_km, 0, zoom_km))
        a.contour(np.flipud(possible_on[z] >= 0), [0.5], colors="white", linewidths=0.5, extent=(0, zoom_km, 0, zoom_km))
        a.set_title(f"Zoom {zoom_km} x {zoom_km} km (cyan box above), {title}\ncyan: burned, white: possibly burned",
                    fontsize=10)
    ax[5].axis("off")
    ax[5].text(0, 0.95, "How to read the zoom\n\n"
               "False colour shows shortwave infrared as red, near infrared as green\n"
               "and red light as blue. Growing crops look green, harvested straw and\n"
               "bare soil pink, char black or dark maroon.\n\n"
               "A field outlined in cyan should look pink or green on the first day\n"
               "and black or dark maroon on the last.", va="top", fontsize=11)
    for a in ax[:2] + ax[3:5]:
        a.set_xlabel("km"); a.set_ylabel("km")
    fig.suptitle("Crop residue burning in the rice belt north-west of Lahore (Sheikhupura, Nankana Sahib, Muridke)", fontsize=13)
    fig.savefig(out, dpi=80)


def main():
    start, end = sys.argv[1], sys.argv[2]
    days = list(burn_data.pass_days(start, end))
    burned_on, possible_on, since, table = detect(days)
    os.makedirs(OUT, exist_ok=True)
    stem = os.path.join(OUT, f"burns_{start}_{end}")
    np.savez_compressed(stem + ".npz", burned_on=burned_on, possible_on=possible_on, since=since, days=np.array(days))
    table.round(3).to_csv(stem + ".csv", index=False)
    figure(days, burned_on, possible_on, table, stem + ".png")
    print(table.round(2).to_string(index=False))
    print(f"burned: {(burned_on >= 0).sum() * HA:,.0f} ha; possibly burned: {(possible_on >= 0).sum() * HA:,.0f} ha")


if __name__ == "__main__":
    main()
