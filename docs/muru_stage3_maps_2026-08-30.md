# Muru redo Stage 3 — the 216-view map suite: COMPLETE

**2026-08-30.** Plan: `validation/muru_redo_plan_2026-08-30.md` §3. Script:
`scripts/muru_redo/stage3_views.py`, run **on nucosmo** (Kev's request;
~2–3 min per halo). Output: `nucosmo:/Volumes/Data2/hestia_muru/stage3_maps/`
(= NFS `/private/nfs/Data2/hestia_muru/stage3_maps/`) — one npz per halo,
~450 MB total, with `PROVENANCE.json` carrying every parameter and the gate
results.

**Per view (36 azimuths × 6 halos), four maps** on the 40°×40°, 0.1° analysis
grid, gcepy orientation, observer at R⊙ = 8.18 kpc in the Stage 1 disk plane,
azimuth from the Stage 1 bar axis:

- `jmap` — correct J = Σ mᵢρᵢ/sᵢ² per pixel / ΔΩ (s ≤ 100 kpc; Stage 2 ρ).
- `stars` / `starsold` — stellar flux Σ mᵢ/sᵢ² /ΔΩ, wind-cleaned; `starsold`
  applies Muru's own old-star cut (a_form < 0.79), their bulge proxy.
- `muru` — their statistic (mass in a 3° flat-(l,b) cone, 15 kpc cut),
  exponent left unapplied (squaring is a monotonic relabeling).

**Why it was cheap enough for nucosmo:** the aperture-KDTree route would have
cost ~18 h there; per-pixel histograms (J, stars) and — because their
aperture is *flat* in (l,b) — an exact disk-kernel convolution (their
statistic) cost seconds per view. The reformulation was **gated, not
assumed**, against the anchored exact implementations on G3.2 view 0:

| gate | median | max | verdict |
|---|---|---|---|
| convolution vs exact `project_muru(metric="flat_lb")` | 0.51% | 2.1% | PASS (rim pixelization; tolerance 1%/5% set above the measurement) |
| histogram J vs aperture `project_annihilation`, 1° smoothing | 0.71% | 15.4% (map edges) | PASS (median gate 5%) |

One gate-side orientation slip (a `.T` carried over from the Julia
comparison, where it was correct) produced a 47% "failure" on first run and
was caught by the gate itself before any production map existed.

## First look — one view (G3.2, α = 0), fractions of ROI total

| map | f(<1°) | f(<3°) | f(<10°) |
|---|---|---|---|
| correct J | 0.077 | 0.312 | 0.734 |
| their (Σm)² | 0.045 | 0.299 | 0.789 |
| old stars | 0.024 | 0.144 | 0.530 |

Two early readings, both to be settled properly in Stage 4 with matched
smoothing, not from this table:

1. **The toy-model contrast (§3 of `muru_projection_analysis.md`: 1–2 orders
   of magnitude in angular scale) is muted on the real simulation** — the
   correct J is more concentrated than their statistic, but mildly. That is
   exactly what §4 predicted: the 0.22 kpc softening (1.5° at 8.2 kpc)
   flattens the cusp precisely where the concentration difference lives, and
   their 3° cone smooths both. The resolution caveat leads every downstream
   statement, as planned.
2. **J vs the old-star map already separates cleanly** (f(<3°) = 0.31 vs
   0.14): even at simulation resolution, the annihilation morphology and the
   old-bulge morphology are not the same object on this view. Whether that
   survives all 216 views and matched smoothing is Stage 4's question.

Azimuthal stability: J's f(<3°) spans 0.279–0.312 over G3.2's 36 views —
the concentration is view-robust at the few-percent level; shape statistics
may vary more.
