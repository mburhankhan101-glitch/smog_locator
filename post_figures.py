"""
Phone-sized slides (1080 x 1350, LinkedIn's portrait format) from the weekly burn map's data, one point
per slide in large type, plus the same slides as one PDF for a LinkedIn document (carousel) post.

    python post_figures.py OUT_DIR

Reads runs/districts/map/burn_map.json (written by burn_map.py); every number comes from there.
"""
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle

DATA = os.path.join("runs", "districts", "map", "burn_map.json")
BG, INK, MUTED, LINE, LAND = "#f3f5f6", "#1b2126", "#56626b", "#c8cfd4", "#dfe4e7"
EMBER, GREY = "#b93d0e", "#8b98a2"
RAMP = ["#f6dcae", "#f1b46c", "#e2783a", "#b93d0e", "#6f1a09"]
BINS = [0.5, 2, 5, 10, 20]                          # ha per km², as on the map page
HEAD, BODY = "Segoe UI", "Segoe UI"                # Segoe UI has a real bold; Bahnschrift does not
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
PITCH_HA = 105 * 68 / 10_000                        # a FIFA-size football pitch
LINKEDIN = "linkedin.com/in/burhan-khan-b4458737b"
GITHUB = "github.com/mburhankhan101-glitch/smog_locator"


def day(s):
    y, m, d = s.split("-")
    return f"{int(d)} {MONTHS[int(m) - 1]}"


def page():
    fig = plt.figure(figsize=(10.8, 13.5), dpi=100, facecolor=BG)
    return fig


def text(fig, x, y, s, size, **kw):
    kw.setdefault("color", INK)
    kw.setdefault("family", BODY)
    return fig.text(x, y, s, fontsize=size, va="top", **kw)


def footer(fig, n, total):
    text(fig, 0.06, 0.045, f"Burhan Khan · {LINKEDIN}", 15, color=MUTED)
    text(fig, 0.94, 0.045, f"{n} / {total}", 15, color=MUTED, ha="right")


def paths(d):
    """SVG path data (M x y L x y ... Z, in km) as a list of (x, y) arrays."""
    out = []
    for part in re.findall(r"M([^Z]*)Z", d):
        pts = np.array([[float(v) for v in p.split()] for p in part.split("L")])
        out.append(pts)
    return out


def map_slide(fig, D, n, total):
    T = D["totals"]
    text(fig, 0.06, 0.96, "Farmland burned\naround Lahore", 58, family=HEAD, weight="bold", linespacing=0.95)
    text(fig, 0.06, 0.838, f"1 Oct to {day(D['latest_image'])} 2026, mapped from Sentinel-2 burn scars",
         20, color=MUTED)
    rows, cols = D["grid"]["rows"], D["grid"]["cols"]
    v = np.zeros(rows * cols)
    for c, w, ha10 in D["cells"]["2026"]:
        v[c] += ha10 / 10
    v = v.reshape(rows, cols)
    v[v < BINS[0]] = np.nan
    wkm, hkm = D["size_km"]
    top, bottom = 0.80, 0.215                         # between the subtitle and the legend
    h = top - bottom
    w = h * 13.5 / 10.8 * wkm / hkm                   # same km per inch across and down
    if w > 0.94:
        w, h = 0.94, 0.94 * 10.8 / 13.5 * hkm / wkm
    ax = fig.add_axes([(1 - w) / 2, top - h, w, h])
    ax.set_facecolor(LAND)
    ax.imshow(v, cmap=ListedColormap(RAMP), norm=BoundaryNorm(BINS + [1e9], 5), extent=(0, wkm, hkm, 0),
              interpolation="nearest")
    for t in D["tehsils"]:
        for p in paths(t["path"]):
            ax.plot(p[:, 0], p[:, 1], color=MUTED if t["reported"] else LINE, lw=0.8,
                    ls="-" if t["reported"] else (0, (3, 3)))
    for d in D["districts"].values():
        for p in paths(d):
            ax.plot(p[:, 0], p[:, 1], color=INK, lw=1.6, alpha=0.8)
    show = {"Lahore City": "LAHORE", "Sheikhupura": "Sheikhupura", "Nankana Sahab": "Nankana Sahib",
            "Pindi Bhattian": "Pindi Bhattian", "Kasur": "Kasur", "Jaranwala": "Jaranwala", "Muridke": "Muridke"}
    for t in D["tehsils"]:
        if t["name"] in show and t["label"]:
            big = t["name"] == "Lahore City"
            ax.text(*t["label"], show[t["name"]], fontsize=19 if big else 15, family=HEAD,
                    weight="bold" if big else "normal", ha="center", va="center", color=INK,
                    bbox=dict(boxstyle="round,pad=0.25", fc=BG, ec="none", alpha=0.85))
    ax.plot([4, 14], [hkm - 4, hkm - 4], color=INK, lw=3)
    ax.text(9, hkm - 5.5, "10 km", ha="center", fontsize=13, color=INK)
    ax.set_xlim(0, wkm), ax.set_ylim(hkm, 0), ax.axis("off")
    labels = ["0.5–2", "2–5", "5–10", "10–20", "20+"]
    y = 0.192
    text(fig, 0.06, y, "Hectares burned per square km:", 17, color=MUTED)
    for i, (c, l) in enumerate(zip(RAMP, labels)):
        x = 0.43 + i * 0.105
        fig.patches.append(Rectangle((x, y - 0.019), 0.022, 0.017, transform=fig.transFigure, fc=c, ec=LINE))
        text(fig, x + 0.028, y, l, 17, color=MUTED)
    ha = T["2026"]["burned_ha"]
    text(fig, 0.06, 0.155, f"{ha:,.0f} ha", 64, family=HEAD, weight="bold", color=EMBER)
    pitches = ha / PITCH_HA
    rounded = round(pitches, -3)
    text(fig, 0.06, 0.088, f"of farmland burned in {len(D['table'])} tehsils, "
         f"{'nearly' if rounded > pitches else 'over'} {rounded:,.0f} football pitches", 22, color=INK)
    footer(fig, n, total)


def missed_slide(fig, D, n, total):
    T = D["totals"]["2026"]
    caught = 1 - T["missed"]
    text(fig, 0.06, 0.955, f"{round(T['missed'] * 100)}%", 190, family=HEAD, weight="bold", color=EMBER)
    text(fig, 0.06, 0.765, "of the burned farmland never\nhad a satellite fire point on it", 40, family=HEAD,
         weight="bold", linespacing=1.05)
    ax = fig.add_axes([0.06, 0.20, 0.88, 0.42])
    filled = round(caught * 100)
    for i in range(100):
        r, c = divmod(i, 10)
        ax.add_patch(Rectangle((c, 9 - r), 0.86, 0.86, fc=EMBER if i < filled else INK, ec="none",
                               alpha=1 if i < filled else 0.82))
    ax.set_xlim(-0.1, 10), ax.set_ylim(-0.1, 10), ax.set_aspect("equal"), ax.axis("off")
    fig.patches.append(Rectangle((0.06, 0.155), 0.03, 0.022, transform=fig.transFigure, fc=EMBER))
    text(fig, 0.105, 0.18, f"Burned, and a fire point saw it: {round(caught * 100)}%", 24)
    fig.patches.append(Rectangle((0.06, 0.112), 0.03, 0.022, transform=fig.transFigure, fc=INK, alpha=0.82))
    text(fig, 0.105, 0.137, f"Burned, and no fire point: {round(T['missed'] * 100)}%", 24)
    text(fig, 0.06, 0.092, f"Each square: 1% of the {T['burned_ha']:,.0f} ha that burned, 1 Oct to {day(T['up_to'])}",
         17, color=MUTED)
    footer(fig, n, total)


def why_slide(fig, D, n, total):
    text(fig, 0.06, 0.955, "Why the fire\npoints miss it", 66, family=HEAD, weight="bold", linespacing=0.95)
    hours = D["viirs_hours"]
    blind = 24 - sum(b - a for a, b in hours)
    ax = fig.add_axes([0.06, 0.62, 0.88, 0.06])
    ax.add_patch(Rectangle((0, 0), 24, 1, fc=LINE))
    for a, b in hours:
        ax.add_patch(Rectangle((a, 0), b - a, 1, fc=EMBER))
    for h in (0, 6, 12, 18, 24):
        ax.text(h, -0.35, f"{h:02d}:00", ha="center", va="top", fontsize=17, color=MUTED)
    ax.set_xlim(0, 24), ax.set_ylim(0, 1), ax.axis("off")
    text(fig, 0.06, 0.735, "When a VIIRS satellite passes (orange)", 22, color=MUTED)
    p = [x for x in D["patches"] if x["season"] == 2026]
    small = [x for x in p if x["patch_pixels"] in ("1", "2-8")]
    caught = sum(x["viirs_compared_ha"] * (1 - x["viirs_missed_share"]) for x in small)
    small_missed = 1 - caught / sum(x["viirs_compared_ha"] for x in small)
    points = [
        (f"Blind about {blind:.0f} hours a day.", "A field that burns in the late afternoon is never seen alight."),
        (f"{round(small_missed * 100)}% of small burns missed.", "Patches under 0.36 ha hide inside a 375 m fire pixel."),
        ("Cloud and smoke block the view.", "The flames are gone; the black field is still there."),
    ]
    y = 0.53
    for head, body in points:
        text(fig, 0.06, y, head, 32, family=HEAD, weight="bold")
        text(fig, 0.06, y - 0.045, body, 23, color=MUTED)
        y -= 0.135
    text(fig, 0.06, 0.135, "A burn scar stays black for days, whatever hour\nthe field burned. Sentinel-2 sees it on the next clear image.",
         22, color=INK, linespacing=1.3)
    footer(fig, n, total)


def tehsil_slide(fig, D, n, total):
    text(fig, 0.06, 0.955, "Where it burned\nmost", 66, family=HEAD, weight="bold", linespacing=0.95)
    text(fig, 0.06, 0.79, f"Share of each tehsil's farmland burned, 1 Oct to {day(D['latest_image'])} 2026", 22, color=MUTED)
    rows = sorted(D["table"], key=lambda r: -r["share_2026"])[:8]
    ax = fig.add_axes([0.38, 0.17, 0.52, 0.58])
    ys = np.arange(len(rows))[::-1]
    ax.barh(ys, [r["share_2026"] * 100 for r in rows], color=[EMBER] + [GREY] * (len(rows) - 1), height=0.62)
    for y, r in zip(ys, rows):
        ax.text(-0.3, y, r["tehsil"], ha="right", va="center", fontsize=24, family=HEAD, color=INK)
        ax.text(r["share_2026"] * 100 + 0.25, y, f"{r['share_2026'] * 100:.1f}%", va="center", fontsize=22,
                family=HEAD, weight="bold", color=INK)
    ax.set_xlim(0, max(r["share_2026"] for r in rows) * 100 * 1.25), ax.axis("off")
    top = rows[0]
    text(fig, 0.06, 0.14, f"{top['tehsil']}: {top['burned_ha_2026']:,.0f} ha burned this season,\n"
         f"and {top['burned_ha_2025_same_date']:,.0f} ha by the same date in 2025.", 21, color=INK, linespacing=1.35)
    text(fig, 0.06, 0.082, "Not yet checked on the ground.", 21, color=MUTED)
    footer(fig, n, total)


def weekly_slide(fig, D, n, total):
    T = D["totals"]
    c = T["2026"]["burned_ha"] / T["2025_same_date"]["burned_ha"] - 1
    text(fig, 0.06, 0.955, "2026 against\n2025", 66, family=HEAD, weight="bold", linespacing=0.95)
    text(fig, 0.06, 0.79, "New burn scars each week (hectares of farmland)", 22, color=MUTED)
    a, b = T["weekly"]["2025"], T["weekly"]["2026"]
    n_weeks = len(b) + 2
    ax = fig.add_axes([0.1, 0.27, 0.84, 0.48])
    x = np.arange(n_weeks)
    ax.bar(x - 0.21, [w["burned_ha"] for w in a[:n_weeks]], 0.4, color=GREY, label="2025")
    ax.bar(x[:len(b)] + 0.21, [w["burned_ha"] for w in b], 0.4, color=EMBER, label="2026")
    ax.bar(x[len(b) - 1] + 0.21, b[-1]["burned_ha"], 0.4, color=BG, ec=EMBER, hatch="///", lw=1.5)
    ax.text(x[len(b) - 1] + 0.21, b[-1]["burned_ha"] + 300, "so far", ha="center", fontsize=16, color=EMBER)
    ax.set_xticks(x, [day(w["start"]) for w in a[:n_weeks]], fontsize=16, color=MUTED, rotation=0)
    ax.tick_params(axis="y", labelsize=16, colors=MUTED)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(LINE)
    ax.grid(axis="y", color=LINE, lw=0.8), ax.set_axisbelow(True), ax.set_facecolor(BG)
    ax.legend(fontsize=20, frameon=False, loc="upper left")
    for lab in ax.get_xticklabels()[1::2]:
        lab.set_visible(False)
    text(fig, 0.06, 0.18, f"By {day(D['latest_image'])}: {T['2026']['burned_ha']:,.0f} ha in 2026, "
         f"{T['2025_same_date']['burned_ha']:,.0f} ha in 2025 ({c:+.0%}).", 24)
    text(fig, 0.06, 0.13, f"Whole 2025 season: {T['2025_season']['burned_ha']:,.0f} ha.", 24, color=MUTED)
    footer(fig, n, total)


def end_slide(fig, D, n, total):
    text(fig, 0.06, 0.955, "See every tehsil,\nweek by week", 66, family=HEAD, weight="bold", linespacing=0.95)
    lines = [("Interactive map", "link in the comments"), ("Code and method", GITHUB),
             ("Data", "Sentinel-2 · ESA WorldCover · NASA FIRMS VIIRS"),
             ("Limits", "Scars are a lower bound; not yet checked on the ground")]
    y = 0.70
    for head, body in lines:
        text(fig, 0.06, y, head, 22, color=MUTED)
        text(fig, 0.06, y - 0.035, body, 24, family=HEAD, weight="bold")
        y -= 0.13
    text(fig, 0.06, 0.16, "Burhan Khan", 40, family=HEAD, weight="bold", color=EMBER)
    text(fig, 0.06, 0.11, LINKEDIN, 24)
    footer(fig, n, total)


def main(out):
    D = json.load(open(DATA, encoding="utf-8"))
    slides = [map_slide, missed_slide, why_slide, tehsil_slide, weekly_slide, end_slide]
    os.makedirs(out, exist_ok=True)
    with PdfPages(os.path.join(out, "smog_locator_carousel.pdf")) as pdf:
        for k, make in enumerate(slides, 1):
            fig = page()
            make(fig, D, k, len(slides))
            fig.savefig(os.path.join(out, f"{k}_{make.__name__.replace('_slide', '')}.png"), facecolor=BG)
            pdf.savefig(fig, facecolor=BG)
            plt.close(fig)
    print(f"{len(slides)} slides in {out}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join("runs", "districts", "post"))
