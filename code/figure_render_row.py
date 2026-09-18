"""One-galaxy, FIVE-panel row figure for the PRL Comment (local render).

Columns are the SAME FIVE as the Supplemental figure (`figure_render.py
--sm`), so the main-text galaxy sits on exactly the same footing as the
other five (Jason Kumar, 2026-09-07). Panel geometry is identical too: the
Supplemental grid uses LEFT/RIGHT margins 0.76/0.04 in at TEXTWIDTH_IN over
five columns, giving 1.256 in square cells -- reproduced here exactly, and
rendered by importing that figure's own `panel()` rather than a copy, so
the two figures cannot drift apart.

PRL counts a two-column figure as 300/(width/height) + 40 word-equivalents.
Adding the fifth column at fixed width SHRINKS the cells, so the row gets
wider in aspect and CHEAPER: aspect 7.08/2.136 = 3.31 costs ~131
word-equivalents, against ~147 for the previous four-column row. It buys
back ~16 words of the 750-word budget rather than spending any.

Galaxy: G1.1 by default -- the halo whose published axis ratios our frame
reconstruction reproduces best (max |ours - theirs| = 0.02 on the DM panel;
the others range to 0.12), and the top row of their Fig. 1, so a reader can
compare panel to panel. NOT because it is strongly barred: its m=2
amplitude is 0.080, second lowest of the six, and the most strongly barred
halo (G2.1, A2 = 0.281) is the WORST match. An earlier note in
analysis/docs/muru_figure_repro_2026-08-30.md claimed a barredness
explanation; it is retracted there -- match quality does not track bar
strength.

Column order and contour-level rule follow figure_render.py exactly:

  1  DM, Ref. [1]             int rho_DM s^2 ds       levels matched to stars
  2  stars, Ref. [1]          int rho_* s^2 ds        levels 70/50/30
  3  "quadratic DM", Ref. [1] [int rho_DM s^2 ds]^2   levels matched to stars
  4  stars                    int rho_* ds            levels 70/50/30
  5  annihilation             J = int rho_DM^2 ds     levels matched to col 4

Columns 1-3 are Ref. [1]'s own three panels in their own conventions,
anchored to their Julia pipeline at 5e-7 (analysis/docs/muru_anchor_2026-08-30.md).
Columns 4-5 are the corrected pair: the s^2 weighting removed from the
stellar map (an old stellar population dilutes as 1/s^2 exactly as
annihilation does -- Jason Kumar, 2026-08-31) and the square taken INSIDE
the line-of-sight integral.

Usage:
    python figure_render_row.py [HALO]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
# `panel` is imported, NOT copied: the main-text row and the Supplemental
# grid must draw identically or the "same footing" claim is cosmetic only.
from figure_render import (NFS, OUT, ST_LEVELS, PUBLISHED, matched_levels,  # noqa: E402
                           panel, TEXTWIDTH_IN, FONT_PT)

HALO = sys.argv[1] if len(sys.argv) > 1 else "G1.1"
STAMP = "2026-09-07"

TITLES = (
    "DM, Ref.~[1]\n" + r"$\propto\!\int\!\rho_{\mathrm{DM}}s^{2}\mathrm{d}s$",
    "stars, Ref.~[1]\n" + r"$\propto\!\int\!\rho_{\star}s^{2}\mathrm{d}s$",
    "``quadratic DM''\n"
    + r"$\bigl[\int\!\rho_{\mathrm{DM}}s^{2}\mathrm{d}s\bigr]^{2}$",
    "stars\n" + r"$\propto\!\int\!\rho_{\star}\mathrm{d}s$",
    "annihilation\n" + r"$J\!=\!\int\!\rho_{\mathrm{DM}}^{2}\mathrm{d}s$",
)
KEYS = ("mu_dm", "stars", "mu_dm_sq", "stars_col", "j_corr")


def main():
    z = np.load(NFS / "figure_maps.npz")
    meta = json.load(open(NFS / "figure_maps_meta.json"))
    OUT.mkdir(parents=True, exist_ok=True)

    mu = z[f"{HALO}_mu_dm"].astype(float)
    st = z[f"{HALO}_mu_stars"].astype(float)          # their cone statistic
    # both corrected maps get THE SAME 3-deg beam their cone statistic
    # carries, so all five panels sit at one resolution
    kern = np.hypot(*np.meshgrid(*(np.arange(-30, 31) * 0.1,) * 2)) <= 3.0
    kern = kern / kern.sum()
    j = fftconvolve(z[f"{HALO}_jmap_muru_geom"], kern, mode="same")
    sc = fftconvolve(z[f"{HALO}_stars_col"], kern, mode="same")

    panels = (mu, st, mu ** 2, sc, j)
    lvls = (matched_levels(mu, st), list(ST_LEVELS),
            matched_levels(mu ** 2, st), list(ST_LEVELS),
            matched_levels(j, sc))
    NCOL = len(panels)

    # Margins in inches, panel size following -- the Supplemental figure's
    # --sm values, so the cells come out the same 1.256 in square. LEFT
    # holds three stacked items: the shared "b [deg]", the "halo + A_2"
    # label, then the y tick numbers.
    LEFT, RIGHT, TOPM, BOTM = 0.76, 0.04, 0.44, 0.44
    W = TEXTWIDTH_IN
    pw = (W - LEFT - RIGHT) / NCOL
    H = pw + TOPM + BOTM
    fig, axes = plt.subplots(1, NCOL, figsize=(W, H), sharey=True,
                             gridspec_kw=dict(wspace=0.0))
    fig.subplots_adjust(left=LEFT / W, right=1 - RIGHT / W,
                        bottom=BOTM / H, top=1 - TOPM / H)

    out = {}
    for k, (ax, m, lv) in enumerate(zip(axes, panels, lvls)):
        out[KEYS[k]] = {
            "levels": lv,
            "ratios": panel(ax, m, lv, k == 0, HALO,
                            meta[HALO]["bar_A2_max"], fs=FONT_PT)}
        ax.set_title(TITLES[k], fontsize=FONT_PT)

    xc = 0.5 * (LEFT + (W - RIGHT)) / W
    yc = 0.5 * (BOTM + (H - TOPM)) / H
    fig.text(xc, 0.004, r"$\ell$ [deg]", ha="center", va="bottom",
             fontsize=FONT_PT)
    fig.text(0.004, yc, r"$b$ [deg]", ha="left", va="center", rotation=90,
             fontsize=FONT_PT)
    for ext in ("pdf", "png"):
        # no bbox_inches="tight" -- keeps the exact printed width
        fig.savefig(OUT / f"muru_fig_row_{HALO}_{STAMP}.{ext}", dpi=200)
    with open(OUT / f"muru_fig_row_{HALO}_{STAMP}.json", "w") as f:
        json.dump({"halo": HALO, "ours": out, "published": PUBLISHED[HALO]},
                  f, indent=1)

    print(f"wrote {OUT}/muru_fig_row_{HALO}_{STAMP}.pdf")
    print(f"  size {W:.2f} x {H:.3f} in, aspect {W/H:.2f}, "
          f"PRL cost ~{300/(W/H) + 40:.0f} word-equivalents, "
          f"cell {pw:.3f} in")
    for key, pub in (("mu_dm", "dm"), ("stars", "stars"), ("mu_dm_sq", "dm2")):
        o = out[key]
        print(f"  {key:<9s} " + "  ".join(
            f"{int(round(lv*100)):>2d}%:{r:.2f}|{pl}%:{pr:.2f}"
            for (lv, r), (pl, pr) in zip(zip(o["levels"], o["ratios"]),
                                         PUBLISHED[HALO][pub])))
    for key in ("stars_col", "j_corr"):
        o = out[key]
        print(f"  {key:<9s} " + "  ".join(f"{int(round(lv*100)):>2d}%:{r:.2f}"
                                          for lv, r in zip(o["levels"],
                                                           o["ratios"])))


if __name__ == "__main__":
    main()
