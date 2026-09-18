"""Render the Muru et al. figure reproduction in THEIR style (local).

Style matched to their plot_methods_densityprojections.jl (read, not
guessed): LINEAR color scale in the :jet1 blue->red map, black contours at
[0.3, 0.5, 0.7, 0.9] of each panel's max, per-contour axis ratios from an
SVD of the REGION above the level (their calc_axis_ratio_svd — grid points
with dens >= level, not the contour curve), printed bottom-right, l axis
reversed (astronomical), a white 2 kpc scale bar of half-width
atan(1 kpc / 8 kpc) = 7.125 deg at b = 16.9 (their accurate value), and a
softening circle — drawn at OUR honest 1.575 deg (0.22 kpc at 8 kpc);
their arc uses 0.533 deg, an apparent h-factor slip in their plotting
code, noted in the results file.

Columns, per Kev 2026-08-30 (their order + the corrected panel):
  1  dark matter        — their cone statistic  ~ int rho_DM s^2 ds
  2  stars (old)        — same statistic on their old-star cut
  3  "quadratic dark matter" — col 1 squared, THEIR fraction levels applied
     to the squared map (as their code does): its q% contour is col 1's
     sqrt(q)% contour — the same curve family relabeled.
  4  corrected annihilation J = int rho_DM^2 ds (same view, 15 kpc,
     3-deg beam to match their aperture's smoothing).

LaTeX fonts (usetex; latex + dvipng present on this machine).
Frames: approximate (stellar shape-tensor major axis) until Muru's AHF
axes arrive — figure_maps.py header.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

#: 🚨 FIGURES ARE BUILT AT FINAL PRINTED SIZE. REVTeX two-column text width
#: is 7.08 in; a figure drawn at 13 in and included with width=\textwidth is
#: scaled by 0.54, shrinking every label to roughly half the caption's font.
#: Drawing at TEXTWIDTH_IN with FONT_PT >= the caption size (PRL captions are
#: ~9 pt) keeps \includegraphics[width=\textwidth] a 1:1 placement.
#: Do NOT pass bbox_inches="tight" for manuscript figures: it crops to the
#: ink and the resulting rescale undoes this.
TEXTWIDTH_IN = 7.08
#: TWO FONT ROLES (Kev, 2026-08-31): axis furniture -- axis labels, tick
#: labels, panel titles -- at 9 pt, matching the caption; in-panel
#: annotations -- the contour level/ratio tables and the scale-bar label --
#: at 6 pt, so they do not crowd the maps.
FONT_PT = 9        # axis labels, tick labels, titles
ANNOT_PT = 6       # contour level/ratio tables, scale-bar label

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath}",
    "font.size": FONT_PT,
    "axes.labelsize": FONT_PT,
    "axes.titlesize": FONT_PT,
    "xtick.labelsize": FONT_PT,
    "ytick.labelsize": FONT_PT,
})

NFS = Path("/private/nfs/Data2/hestia_muru")
OUT = Path(__file__).resolve().parents[1] / "figures"
HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
#: 🚨 THE PARTICLE FILES ARE MISLABELLED WITHIN EACH PAIR relative to Ref. [1]
#: (found 2026-09-16; analysis/docs/muru_ahf_frames_2026-09-15.md Finding 7).
#: Every key in this codebase -- HALOS, figure_maps.npz, muru_frames.json,
#: stage3/4 outputs -- is the PARTICLE-FILE label. The paper's Table I label
#: for the same halo is the other member of the pair, anchored two ways:
#: M200c + R200c of the AHF profiles match Table I by name, and each profile
#: matches the OTHER particle file on M(<r) to 0.2%. Rows are drawn in PAPER
#: order and labelled with PAPER names; PUBLISHED is looked up by paper name.
#: Do not "fix" this by renaming the data keys -- too many products carry
#: them; map at the presentation layer only.
PAPER_LABEL = {"G1.1": "G1.2", "G1.2": "G1.1", "G2.1": "G2.2",
               "G2.2": "G2.1", "G3.1": "G3.2", "G3.2": "G3.1"}
ST_LEVELS = (0.7, 0.5, 0.3)     # published App. A: stellar contours fixed;
#                                 DM/DM2 levels derived by b=0 l-intercept
#                                 matching against the stellar contours
X = np.linspace(-19.95, 19.95, 400)

#: Published printed corner-table values, READ FROM THE PRL PDF (rule 1a's
#: approved channel — these are typeset numbers, not measurements off
#: curves): Figs. 1-2 of PRL 135, 161005 (2025). Format:
#: {halo: {"dm"/"stars"/"dm2": [(level_pct, ratio), ...]}}. Note their DM2
#: levels are their DM levels squared and carry IDENTICAL ratios — their
#: own table demonstrates the monotonic relabeling.
PUBLISHED = {
    "G1.1": {"dm": [(69, .82), (52, .76), (37, .68)],
             "stars": [(70, .63), (50, .58), (30, .48)],
             "dm2": [(47, .82), (27, .76), (13, .68)]},
    "G2.1": {"dm": [(74, .79), (60, .72), (43, .64)],
             "stars": [(70, .68), (50, .62), (30, .52)],
             "dm2": [(55, .79), (36, .72), (18, .64)]},
    "G3.1": {"dm": [(70, .74), (51, .68), (29, .73)],
             "stars": [(70, .54), (50, .44), (30, .35)],
             "dm2": [(49, .74), (26, .68), (9, .73)]},
    "G1.2": {"dm": [(75, .82), (58, .78), (39, .71)],
             "stars": [(70, .77), (50, .72), (30, .62)],
             "dm2": [(56, .82), (34, .78), (16, .71)]},
    "G2.2": {"dm": [(73, .81), (58, .77), (40, .73)],
             "stars": [(70, .71), (50, .63), (30, .55)],
             "dm2": [(53, .81), (34, .77), (16, .73)]},
    "G3.2": {"dm": [(76, .76), (58, .69), (41, .65)],
             "stars": [(70, .66), (50, .59), (30, .50)],
             "dm2": [(57, .76), (34, .69), (17, .65)]},
}
SCALEBAR_HALF = np.degrees(np.arctan(1.0 / 8.0))  # 2 kpc at 8 kpc = 7.125 deg
SOFT_DEG = np.degrees(np.arctan(0.22 / 8.0))      # honest softening: 1.575 deg


def region_axis_ratio(m, level):
    """Their calc_axis_ratio_svd: SVD of the centered coords of grid points
    with dens >= level; min/max singular value."""
    ib, il = np.where(m >= level * m.max())
    if len(ib) < 8:
        return np.nan
    c = np.column_stack([X[il], X[ib]])
    c = c - c.mean(axis=0)
    s = np.linalg.svd(c, compute_uv=False)
    return float(s[1] / s[0])


def b0_profile(m):
    """The map along b = 0 (GC on the row 199/200 boundary)."""
    return 0.5 * (m[199] + m[200])


def matched_levels(m, st):
    """Their published rule: contour levels of map `m` chosen so its contours
    cross b = 0 at the same +l as the stellar 70/50/30 contours."""
    ps, pm = b0_profile(st / st.max()), b0_profile(m / m.max())
    out = []
    for L in ST_LEVELS:
        idx = np.where(ps >= L)[0].max()          # rightmost (+l) crossing
        out.append(float(pm[idx]))
    return out


def panel(ax, m, levels, first_col, halo, a2, fs=FONT_PT,
          fa=ANNOT_PT):
    m = m / m.max()
    # aspect="auto": the axes boxes are laid out square by construction
    # below, so the image fills each cell exactly and the panels abut
    ax.imshow(m, origin="lower", extent=(-20, 20, -20, 20), aspect="auto",
              vmin=0.0, vmax=1.0, cmap="jet", interpolation="nearest")
    ax.contour(X, X, m, levels=sorted(levels), colors="k", linewidths=0.9)
    ratios = [region_axis_ratio(m, lv) for lv in levels]
    # header dropped: at >=9 pt it spanned the panel and overlapped the
    # contours; the caption carries "level; axis ratio" instead
    txt = r"\begin{tabular}{@{}r@{\ \ }r@{}}"
    for lv, r in zip(levels, ratios):
        rr = "--" if np.isnan(r) else f"{r:.2f}"
        txt += rf"{int(round(lv*100))}\% & {rr}\\"
    txt += r"\end{tabular}"
    ax.text(0.97, 0.03, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=fa, color="w")
    # accurate 2 kpc scale bar (their placement)
    ax.plot([-SCALEBAR_HALF, SCALEBAR_HALF], [16.4, 16.4], color="w", lw=1.2)
    ax.text(0, 17.0, r"2\,kpc", ha="center", va="bottom", fontsize=fa,
            color="w")
    # softening circle, honest subtense
    ax.add_patch(plt.Circle((-16.5, 16.8), SOFT_DEG, fill=False,
                            color="w", lw=0.8))
    if first_col:
        # 🚨 A2 lives in the ROW LABEL, outside the panel. Inside it collided
        # with the scale bar and with the level/ratio table once five columns
        # shrank the panels. It is a SINGLE line: the earlier two-line form
        # ("halo + A_2" over "b [deg]") let the A_2 subscript descend into the
        # b below it. The axis names are drawn once for the whole grid.
        ax.set_ylabel(rf"{PAPER_LABEL[halo]}\ \ $A_2={a2:.2f}$", fontsize=fs, labelpad=2)
    # ticks off the panel edges: with zero spacing, edge labels from adjacent
    # panels would overlap
    ax.set_xticks([10, 0, -10])
    ax.set_yticks([-10, 0, 10])
    ax.tick_params(labelsize=fs)
    ax.invert_xaxis()                             # their xreversed
    ax.set_xlim(20, -20)
    return ratios


def main():
    z = np.load(NFS / "figure_maps.npz")
    meta = json.load(open(NFS / "figure_maps_meta.json"))
    OUT.mkdir(parents=True, exist_ok=True)

    from scipy.signal import fftconvolve
    kern = np.hypot(*np.meshgrid(*(np.arange(-30, 31) * 0.1,) * 2)) <= 3.0
    kern = kern / kern.sum()

    # Five columns: the three panels of Ref. [1] as they plot them, then the
    # flux-weighted pair. The stellar column density int rho_* ds -- not their
    # s^2-weighted cone mass -- is the like-for-like comparison for the
    # annihilation map, since a stellar source dilutes as 1/s^2 too
    # (Jason Kumar, 2026-08-31).
    titles = (
        "DM, Ref.~[1]\n"
        + r"$\propto\!\int\!\rho_{\mathrm{DM}}s^{2}\mathrm{d}s$",
        "stars, Ref.~[1]\n"
        + r"$\propto\!\int\!\rho_{\star}s^{2}\mathrm{d}s$",
        "``quadratic DM''\n"
        + r"$\bigl[\int\!\rho_{\mathrm{DM}}s^{2}\mathrm{d}s\bigr]^{2}$",
        "stars\n" + r"$\propto\!\int\!\rho_{\star}\mathrm{d}s$",
        "annihilation\n" + r"$J\!=\!\int\!\rho_{\mathrm{DM}}^{2}\mathrm{d}s$",
    )
    NCOL = len(titles)
    # --sm omits the galaxy carried in the main text (default G1.1), so the
    # Supplemental figure shows the REMAINING galaxies rather than repeating it.
    #
    # --all is the COMBINED figure (Kev, 2026-09-08): all six halos x five
    # columns in ONE figure. PRL's 750-word budget cannot carry a main-text
    # figure as well, so the main-text row was dropped and the galaxy it
    # carried (G1.1) rejoins the SM grid rather than being repeated across two
    # figures. Same SM layout, same cell size (1.256 in), six rows instead of
    # five -- so it stays on the identical footing panel() enforces.
    every = "--all" in sys.argv
    sm = "--sm" in sys.argv or every
    skip = None
    if "--skip" in sys.argv:
        skip = sys.argv[sys.argv.index("--skip") + 1]
    elif sm and not every:
        skip = "G1.1"
    # Row order = Ref. [1]'s figure order (Kev, 2026-09-16): their Fig. 1 is
    # G1.1, G2.1, G3.1 and their Fig. 2 is G1.2, G2.2, G3.2, so the reader can
    # run down our rows against their two figures in sequence.
    PAPER_ORDER = ("G1.1", "G2.1", "G3.1", "G1.2", "G2.2", "G3.2")
    halos = tuple(sorted((h for h in HALOS if h != skip),
                         key=lambda h: PAPER_ORDER.index(PAPER_LABEL[h])))

    # CONTIGUOUS GRID at final printed size. Margins are given in inches and
    # the panel size follows, so every cell is exactly square and the maps
    # abut with no gutters; the axis names are drawn once, at the margins.
    fs = FONT_PT
    # LEFT holds three stacked items: the shared "b [deg]" at the far edge,
    # then the per-row "halo + A_2" label, then the y tick numbers
    LEFT, RIGHT, TOPM, BOTM = (0.76, 0.04, 0.44, 0.44) if sm else \
                              (0.80, 0.04, 0.46, 0.44)
    W = TEXTWIDTH_IN if sm else 2.9 * NCOL
    #: 🚨 SIX ROWS DO NOT FIT AT FULL TEXT WIDTH. At W = 7.08 in the grid is
    #: 606 pt tall against REVTeX's 672 pt \textheight, and float separation
    #: eats the rest: measured, only ~100 caption words fit before LaTeX
    #: reports "Float too large for page". At 0.85 x \textwidth the full
    #: 257-word caption fits with room to spare (measured by compiling, see
    #: analysis/docs/fig_sm_6halo_caption_2026-09-08.md).
    #:
    #: The figure is still drawn at FINAL PRINTED SIZE -- 6.02 in wide, to be
    #: included as \includegraphics[width=0.85\textwidth], which is a 1:1
    #: placement. Do NOT render at 7.08 and let \includegraphics shrink it:
    #: that rescales every label below the caption's font size, which is the
    #: whole point of TEXTWIDTH_IN above.
    if every:
        W = 0.85 * TEXTWIDTH_IN
    if "--width" in sys.argv:
        W = float(sys.argv[sys.argv.index("--width") + 1])
    pw = (W - LEFT - RIGHT) / NCOL
    H = pw * len(halos) + TOPM + BOTM
    fig, axes = plt.subplots(len(halos), NCOL, figsize=(W, H),
                             sharex=True, sharey=True,
                             gridspec_kw=dict(wspace=0.0, hspace=0.0))
    fig.subplots_adjust(left=LEFT / W, right=1 - RIGHT / W,
                        bottom=BOTM / H, top=1 - TOPM / H)
    ratios = {}
    for i, h in enumerate(halos):
        mu = z[f"{h}_mu_dm"].astype(float)
        st = z[f"{h}_mu_stars"].astype(float)
        j = fftconvolve(z[f"{h}_jmap_muru_geom"], kern, mode="same")
        sc = fftconvolve(z[f"{h}_stars_col"], kern, mode="same")
        panels = (mu, st, mu ** 2, sc, j)
        lvls = (matched_levels(mu, st), list(ST_LEVELS),
                matched_levels(mu ** 2, st), list(ST_LEVELS),
                matched_levels(j, sc))
        ratios[h] = {}
        for k, (m, lv) in enumerate(zip(panels, lvls)):
            ax = axes[i, k]
            key = ("mu_dm", "stars", "mu_dm_sq", "stars_col", "j_corr")[k]
            ratios[h][key] = {"levels": lv,
                              "ratios": panel(ax, m, lv, k == 0, h,
                                              meta[h]["bar_A2_max"], fs=fs)}
            if i == 0:
                ax.set_title(titles[k], fontsize=fs)

    # --sm: Supplemental Material version -- NO suptitle (the SM text carries
    # the description; a title inside the figure would duplicate it).
    xc = 0.5 * (LEFT + (W - RIGHT)) / W
    yc = 0.5 * (BOTM + (H - TOPM)) / H
    fig.text(xc, 0.004, r"$\ell$ [deg]", ha="center", va="bottom", fontsize=fs)
    fig.text(0.004, yc, r"$b$ [deg]", ha="left", va="center", rotation=90,
             fontsize=fs)
    if sm:
        # 2026-09-15: AHF frames (Muru's profile files). The 2026-09-08 stem
        # is the pre-AHF reconstruction and is kept for comparison.
        stem = ("muru_fig_sm_6halo_2026-09-16" if every
                else "muru_fig_sm_5halo_2026-08-31")
        # --tag <t>: suffix the stem, so alternative frame choices (e.g. the
        # option-(b) best-azimuth frames, tag "bestaz") do not overwrite the
        # observer-on-Ea render of the same date.
        if "--tag" in sys.argv:
            stem += "_" + sys.argv[sys.argv.index("--tag") + 1]
    else:
        fig.suptitle(
            r"Muru et al.\ (PRL 135, 161005) Figs.~1--2 reproduction "
            r"(their App.~A frames and contour-level rule) $+$ corrected "
            r"annihilation panels. "
            "\n"
            r"Squaring after projection relabels levels only: the $q\%$ contour "
            r"of $\bigl[\int\rho\,s^{2}\mathrm{d}s\bigr]^{2}$ is the "
            r"$\sqrt{q}\%$ contour of $\int\rho\,s^{2}\mathrm{d}s$, whereas the "
            r"annihilation signal is $\int\rho^{2}\,\mathrm{d}s \neq "
            r"\bigl[\int\rho\,s^{2}\,\mathrm{d}s\bigr]^{2}$. "
            r"Circle: softening (1.6$^\circ$); bar: 2 kpc.",
            fontsize=9)
        stem = "muru_figure_repro_2026-08-30"
    for ext in ("pdf", "png"):
        # no bbox_inches="tight": the SM figure must keep its exact printed
        # width so \includegraphics[width=\textwidth] does not rescale fonts
        fig.savefig(OUT / f"{stem}.{ext}", dpi=170)
    for h in ratios:
        ratios[h]["paper_label"] = PAPER_LABEL[h]
    with open(OUT / f"{stem}_ratios.json", "w") as f:
        json.dump(ratios, f, indent=1)
    print("wrote", OUT / f"{stem}.pdf")
    # side-by-side against the PUBLISHED printed values (levels%, ratio)
    print("\n=== ours vs published (level%: ratio | their level%: ratio) ===")
    print("    rows: particle-file key -> PAPER label (published values looked up by paper label)")
    for h in halos:
        for key, pub in (("mu_dm", "dm"), ("stars", "stars"),
                         ("mu_dm_sq", "dm2")):
            ours = ratios[h][key]
            cells = "  ".join(
                f"{int(round(lv*100)):>2d}%:{r:.2f}|{pl}%:{pr:.2f}"
                for (lv, r), (pl, pr) in zip(
                    zip(ours["levels"], ours["ratios"]), PUBLISHED[PAPER_LABEL[h]][pub]))
            print(f"  file {h} = paper {PAPER_LABEL[h]} {key:<9s} {cells}")


if __name__ == "__main__":
    main()
