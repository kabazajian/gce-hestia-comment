# Muru redo Stage 5 — systematics budget: no conclusion flips; one headline demoted

**2026-08-30.** Plan §5 of `validation/muru_redo_plan_2026-08-30.md`.
Scripts: `stage5_systematics.py` (variant maps + statistics on nucosmo,
7.2 min), `stage5_summary.py` (deltas), `stage4b_continuum.py` (now takes
FWHM as an argument for the smoothing axis). Variants on G3.2's 36 views
(the anchored halo) against an identical-code-path baseline; center axis on
G3.1 (the worst Stage 1 shift). All deltas are median (max) |Δ| over views
at FWHM 2°.

## The budget

| axis | corr_gnfw2 | q30 | c4_30 | conc_3_10 | pref. flips |
|---|---|---|---|---|---|
| kNN k = 16 | 0.004 (0.004) | 0.004 (0.010) | 0.001 (0.004) | 0.004 (0.004) | 0/36 |
| kNN k = 64 | 0.001 (0.002) | 0.001 (0.008) | 0.001 (0.004) | 0.002 (0.003) | 0/36 |
| s_max 50 vs 100 kpc | **0.0000** | 0.0000 | 0.0000 | 0.0000 | 0/36 |
| R⊙ 8.0 vs 8.18 | 0.003 (0.004) | 0.002 (0.007) | 0.000 (0.002) | 0.004 (0.004) | 0/36 |
| half-sample a | 0.011 (0.016) | **0.031 (0.064)** | 0.008 (0.015) | 0.008 (0.012) | 0/36 |
| half-sample b | 0.003 (0.010) | **0.033 (0.048)** | 0.006 (0.016) | 0.002 (0.006) | 0/36 |
| center: AHF vs refined (G3.1) | **0.118 (0.192)** | 0.003 (0.010) | 0.011 (0.023) | 0.036 (0.067) | 0/36 |

**Zero preference flips on any axis.** Readings:

1. **k = 32 is confirmed** at the statistics level (≤ 0.005 on everything).
2. **s_max is a non-axis**: 50 vs 100 kpc changes nothing to 4 decimals —
   the ρ² weighting extinguishes the line of sight long before 50 kpc.
3. **R⊙ 8.0 vs 8.18 is negligible** (≤ 0.004).
4. **Particle noise sets the floor**: half-sampling moves q30 by ±0.03
   (max 0.06) and c4 by ±0.01 (max 0.016). Against it: the DM-vs-star q30
   separation (0.76 vs 0.47) is ~5–10 noise units — safe. The tightest
   margin in the whole analysis is the single most-boxy view's c4 (−0.008)
   vs f98's (−0.0242): a 0.016 gap, equal to the *max* observed
   half-sample shift. The "no view as boxy as F98" statement is a
   distribution-level conclusion (median +0.029, ~20 noise units away);
   for the one extreme view it holds at ~1 noise unit and is stated that
   way, not stronger.
5. **Centering is the one large correlation systematic**: G3.1's 0.27 kpc
   AHF-vs-refined offset (≈ 1.9° at 8.18 kpc) moves template correlations
   by median 0.12 — while barely touching the shape statistics (q30 0.003,
   c4 0.011). Correlation statistics need the refined centers, which is
   what every production stage uses; shape statistics are robust to even
   the worst center error we measured.

## The smoothing axis — the one place a headline moved

Re-running the full Stage 4b analysis at FWHM 3° (both smoothings were in
`stage4_stats.json` all along): **the Pearson preference counts are not
smoothing-stable.** Combos 3, 4, 5, 7 (NB count-fractions 0.27–0.30) fall
from 216/216 gNFW²-wins to 115/127/125/67 — heavier smoothing pushes the
comparison further onto the concentration axis, where mid-NB combos
interpolate to the simulated maps' concentration (the combo_10 mechanism,
now seen mid-family). Everything shape-based is unchanged at 3°: f98's c4
rank stays 0.000 (no view as boxy), the DM/star separation stands, the
continuum still fails at leave-one-halo-out (71/216) and
leave-one-realization-out (70/216; G3 pair 0/72), mirror invariance holds.

**Consequence, promoted to a rule for the write-up:** Pearson-preference
counts are secondary evidence (concentration-dominated and
smoothing-fragile); the shape statistics (q, c4) and the continuum result
carry the conclusions. Stage 4's record already flagged this as a caveat —
Stage 5 upgrades it to a measured requirement.

## Standing restrictions (unchanged)

Scales ≥ FWHM 2° ≥ 1.26× the softening subtense; the inner ~1 kpc is
unresolved; frames are ours (AHF axes pending — their arrival tightens the
figure reproduction, not these statistics, which marginalize over azimuth).

**Stage 5 verdict: every Stage 4 conclusion survives every axis tested,
with the Pearson counts demoted as above. The redo's results are ready for
Stage 6 routing (Comment 4.3 / paper §12.D / separate paper — Kev's call).**
