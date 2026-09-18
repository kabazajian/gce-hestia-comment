# Muru redo Stage 4 — morphology statistics, typicality, and the continuum test

**2026-08-30.** Plan §4/§4b of `validation/muru_redo_plan_2026-08-30.md`.
Scripts: `stage4_templates.py` (template construction, Studio — needs repo
data), `stage4_stats.py` (statistics over all 216 views + templates, run on
nucosmo, 4.4 min), `stage4b_continuum.py` (hold-out + typicality, local).
Artifacts beside the data: `stage4_templates.npz` + meta,
`stage4_stats.json`. All statements below are counting statements over the
216 views (6 halos × 36 azimuths) at FWHM 2° smoothing unless noted.

## The comparison templates

Top-10 (bulge+NB, diffuse) combos ranked by ln_H **within the complete-data
cell** `rfree10_nomask`, no plane mask. ⚠ A pooled ranking across masking
cells was tried first and returned ten `rfree0_mask90` rows — cells that
delete the most pixels trivially have the highest ln_H, so cross-cell ln_H
ranking measures data volume, not model quality; recorded in the script as a
guard comment. The top cells are dominated by the Galp21 diffuse families;
all three bulges appear; NB count-fractions 0.13–0.54. The A42 λ-composite
is still the eventual referent (unranked as of 2026-08-23); these stand in.
Pure-morphology maps (Gordon F98/Coleman20, Pohl NB, analytic gNFW² at
R⊙ = 8.18) with the mirror-orientation check passing (0.999 vs 0.985); cao13
combos use gcepy's baked template (stated caveat).

## Result 1 — the stellar side reproduces the bulge; the DM side does not

- **Old-star maps** (their own a_form < 0.79 cut) correlate with the real
  bulge templates at **median r = 0.960 (coleman20) / 0.944 (f98)** against
  0.723 for gNFW². HESTIA's old stellar populations *are* bulge-shaped —
  the pipeline reproduces the stellar-morphology expectation end to end.
- **Corrected J-maps** correlate best with the spherical/concentrated side:
  vs 9 of the 10 combos, **gNFW² wins in 216 of 216 views (6 of 6 halos)**
  (the 10th, combo_10, is the concentration coincidence below).

## Result 2 — shape statistics: the corrected annihilation maps are not boxy

| | q30 (axis ratio) | c4_30 (boxiness; <0 = boxy) | conc f(<3°)/f(<10°) |
|---|---|---|---|
| 216 J views | median 0.76 [0.51, 0.91] | median **+0.029** [−0.008, +0.099] | 0.43 [0.34, 0.50] |
| 216 old-star views | median 0.47 [0.17, 1.00] | median +0.114 [−0.08, +0.46] | — |
| f98 template | 0.715 | **−0.0242** | 0.267 |
| coleman20 | 0.615 | +0.0259 | 0.292 |
| gNFW² | 1.000 | 0.000 | 0.491 |

- **No simulated corrected annihilation map is as boxy as F98**: every one
  of 216 views has c4 above f98's −0.024 (percentile rank 0.000). The
  corrected maps are mildly flattened (halo medians q30 0.61–0.84) and
  slightly *disky*, not boxy — the §5 prediction of
  `muru_projection_analysis.md` ("rounder, not boxier"), now measured on
  their own halos.
- **DM and stars separate within the same view**: J q30 median 0.76 vs
  stars 0.47; c4 +0.03 vs +0.11. With the correct functional on each side
  (ρ² vs ρ), the two morphologies are different objects — their
  "indistinguishable" conclusion does not survive the corrected estimator,
  at the scales the simulation resolves.

## Result 3 — the continuum test FAILS at the load-bearing rung (Kev's 4b)

Gaussian fit to the 5-dim statistic vector; held-out coverage of the 90%
Mahalanobis region:

| rung | inside 90% | expected |
|---|---|---|
| a. random 10% ×20 | 397 / 440 | ~396 ✓ |
| b. 40° arcs (1/halo) | 19 / 24 | ~22 (≈) |
| c. **leave-one-halo-out** | **78 / 216** | ~194 ✗ (G2.2, G3.2: 0/36) |
| d. leave-one-realization-out | 75 / 216 | ~194 ✗ (G3 pair: 0/72) |

Azimuthal smoothness holds (rungs a–b); **the 216 views are six
halo-specific morphology clusters, not draws from one continuum**. Effective
sample size is 3 Local Group realizations. Consequences, per the plan: no
fitted-distribution typicality is quoted (per-statistic counting ranks
only), and the failure is itself a finding — the simulated GCE morphology
is halo-history-specific, so "constrained simulations predict the GCE DM
morphology" is not supportable as a universal statement from six halos.

## Two honest caveats that must travel with all of this

1. **Pearson on linear maps is concentration-dominated.** combo_10
   (17% NB) has conc 0.429 — matching the view median 0.427 almost exactly —
   and beats gNFW² in Pearson in 216/216 views *for both J and stellar
   maps*. That is concentration coincidence, not bar morphology (its shape
   stats sit with the other combos). Similarly J "prefers f98 over gNFW²"
   in 89/216 views by Pearson while being less boxy than f98 in 216/216 by
   c4. Shape statistics carry the morphology conclusions; correlations are
   reported with this caveat attached.
2. **Resolution.** The 0.22 kpc softening (1.58° at 8.2 kpc) flattens the
   inner cusp; all statistics are at FWHM ≥ 2°, and the muted
   concentration contrast between correct-J and their statistic (Stage 3
   first look) is partly the softening. The inner ~1 kpc is unresolved and
   nothing here claims otherwise.

Mirror check: all statistics mirror-invariant to ≤ 0.004 (the F98 bar tilt
contributes the 0.004 on corr_combo) — the 36 azimuths per halo are ~18
independent views, stated for the counting denominators.

## Where this feeds

Comment item 4.3 (protective calculation): **their physical conclusion does
not survive the corrected estimator on their own halos** — corrected maps
are non-boxy and shape-separated from the stellar maps, while their own
old-star cut *does* reproduce the bulge templates. Paper §12.D gets the
one-paragraph version. Stage 5 (systematics: smoothing scale, k, s_max,
center, half-sample) remains before any of this is quoted outside the repo.
