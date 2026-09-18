# Muru et al. Figs. 1–2: reproduced (approximate frames) + corrected panels

**2026-08-30.** Kev's ask: reproduce their figures with the incorrect
projections, exactly or approximately, plus corrected versions. Scripts:
`scripts/muru_redo/figure_maps.py` (map generation, nucosmo, 80 s),
`figure_render.py` (local matplotlib). Figure:
`results/pdf/muru_figure_repro_2026-08-30.{pdf,png}`; per-contour ratios in
`muru_figure_repro_ratios_2026-08-30.json`.

## UPDATE 2026-08-30 (published PRL in hand): frames now per their App. A

The published version (`docs/Muru_PRL_2025.pdf`, PRL 135, 161005) states
the frame construction the arXiv quotes lacked: inertia tensor from ALL
particles ≤ 10 h⁻¹ kpc of the center, observer on its major axis — and the
DM contour-level rule (levels chosen to cross b = 0 at the same ℓ as the
stellar 70/50/30 contours). Both are now implemented (gas included in the
tensor, AHF center, their level rule), giving the first **numerical**
comparison against their printed corner tables:

- **G1.1 reproduces their printed ratios to ≤ 0.02** (DM 0.80/0.76/0.69 vs
  their 0.82/0.76/0.68) with levels within ~5 percentage points.
- The other halos differ by up to ~0.2, **dominated by the stellar
  panels** — since the statistic is anchored, the residual is the tensor
  axes: AHF's exact shape algorithm (weighting/iteration) vs our plain
  mass-weighted second-moment tensor. Muru's AHF profile files (requested)
  settle it.

  🚨 **RETRACTED 2026-08-31 — the barredness explanation was wrong.** This
  section originally read "*Rounder halos have less-determined axes, which
  is why the strongly-barred G1.1 agrees and others drift*." G1.1's m=2
  amplitude is **0.080, the second lowest of the six**, and the most
  strongly barred halo (G2.1, A2 = 0.281) is the **worst** match (0.12).
  Match quality does not track bar strength:

  | halo | bar A2 | max \|ours − theirs\| (DM panel) |
  |---|---|---|
  | G1.1 | 0.080 | **0.02** |
  | G1.2 | 0.149 | 0.06 |
  | G2.1 | **0.281** | **0.12** |
  | G2.2 | 0.061 | 0.11 |
  | G3.1 | 0.079 | 0.11 |
  | G3.2 | 0.194 | 0.07 |

  The correct statement is that **one halo (G1.1) matches to 0.02 and the
  others differ by up to 0.12, attributable to the inertia-tensor
  construction** — with no mechanism claimed. ⚠ The retracted phrasing had
  already propagated into the Comment manuscript ("*for the most strongly
  barred galaxy the printed axis ratios are reproduced to 0.02*"); flagged
  to Kev 2026-08-31 for removal.
- Our derived quadratic-DM levels come out ≈ the DM levels squared with
  identical ratios — as in **their own published corner tables** (e.g.
  G1.1: 69/52/37% → 0.82/0.76/0.68 and 47/27/13% → 0.82/0.76/0.68), which
  therefore demonstrate the monotonic relabeling with no recomputation.
- The published figures carry **no softening circle**; the 0.533° arc (the
  apparent h-slip noted below) is in an unused single-panel code path
  only. Our figure keeps a circle at the honest 1.575° as our own
  addition, so the note below is about their code, not their figures.

Section below kept as the pre-PRL record:

## Fidelity: approximate now, exact when Muru's AHF files arrive

Everything except the frame is **their pipeline exactly** — the cone-mass
statistic, 3° flat-(l,b) aperture, 15 kpc cut, max normalization, and
per-contour SVD axis ratios, all anchored to their Julia code at 5×10⁻⁷
(`validation/muru_anchor_2026-08-30.md`). The one approximation: their
observer sits along the **AHF-profile inertia-tensor major axis** (private
loader; files requested), ours along the in-plane stellar shape-tensor
major axis — which differs from the bar axis by 5°–79° depending on halo,
so the per-panel comparison will sharpen when their axes arrive. Their
printed ratios are graphical values and are compared **qualitatively only**
(rule 1b): our columns reproduce the same extended, moderately flattened
morphology family with ratios declining outward, as their figures show.

**Should our per-contour values match theirs? Not yet — and the size of the
expected mismatch is measured, not assumed** (Kev's question, 2026-08-30).
Viewing azimuth alone moves their statistic's region-SVD ratio by 0.2–0.33
at the 30% level *within a single halo* (measured on the Stage 3 36-azimuth
`muru` maps: G3.2 spans 0.32–0.52, G2.1 spans 0.38–0.71; the 90% level
spans similarly). Since our stand-in axis differs from their AHF axis by an
unknown angle (and from the bar axis by 5°–79°), per-contour differences of
this size are exactly what the frame mismatch predicts. Every other
ingredient is anchored to their code at 5×10⁻⁷, so once their AHF axes
arrive the values should agree to their Float32 — and a residual mismatch
at that point would be a genuine flag to chase, not noise.

## Style (matched to their plotting code, 2026-08-30 restyle per Kev)

Read from their `plot_methods_densityprojections.jl`, not guessed: LINEAR
color scale in `:jet1` (blue→red), black contours at [0.3, 0.5, 0.7, 0.9]
of each panel's max, axis ratios from an SVD of the **region** above the
level (their `calc_axis_ratio_svd` — grid points, not contour vertices),
printed bottom-right as their "Cont. level / Ratio" table, ℓ axis reversed,
white 2 kpc scale bar of half-width atan(1/8) = 7.125° at b = 16.9°.
LaTeX fonts and equations throughout; column titles carry the integrals.

⚠ One discrepancy found while matching: their softening arc is drawn at
0.533° where 0.22 kpc at 8 kpc subtends **1.575°** — an apparent h-factor
slip in their plotting code (their comment converts 220 pc → 148.94 pc/h
and then halves it). Our circle is drawn at the honest 1.575°. Display
cosmetics only — it does not touch their analysis — but worth one line if
the resolution point appears in the Comment.

## What the figure demonstrates

Four columns per halo in their order plus the corrected panel — their
linear DM, their old-star panel, their "quadratic DM" (col 1 squared, their
fraction levels applied to the squared map exactly as their code does), and
the corrected J (same view, same 15 kpc truncation, matched 3° beam):

1. **The "quadratic DM" column is the linear DM column relabeled.** With
   their fraction levels applied to each panel's own max (their code's
   convention), the quadratic panel's q% contour is the linear panel's
   √q% contour — verified in the printed ratios: e.g. G1.1's quadratic 30%
   ratio (0.75) sits exactly between the linear panel's 50% (0.73) and 70%
   (0.79), as √0.30 = 55% demands, in every halo. The squaring adds no
   morphological information — the Comment's monotonicity point on their
   own galaxies.
2. **The corrected J is a different object**: compact and near-round —
   ratios 0.71–0.95 across their four levels, visibly a small round source
   where their DM statistic fills the 40°×40° panel (30%-level ratios
   0.57–0.67, extended to |ℓ| ≳ 15°).
3. **The stellar panels are unmistakably flattened disks** (30%-level
   ratios 0.27–0.58), so the "DM and stellar morphologies are
   indistinguishable" comparison fails visually as well as statistically
   once the DM side is computed correctly.

Physical detail: G2.2's corrected panel shows a secondary compact source
low-right — a DM subhalo inside the 15 kpc truncation along that sight
line (ρ² weighting makes subhalos visible where their mass-sum statistic
buries them). Real structure, and its presence slightly distorts that
panel's 2% ratio (0.96).

## Status

This is the Comment's candidate figure (item §3.5 of
`validation/muru_comment_prep.md` favored the spherical-halo table; this
six-halo version is stronger — same demonstration on their own galaxies).
When the AHF files arrive, swap `figure_maps.major_axis()` for their
delivered axes and rerun both scripts — everything downstream is unchanged.
