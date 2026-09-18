# The gNFW² concentration C = 0.49 — what it is, and a 3.4% quadrature bias

**2026-09-08.** Two findings, both about the single number `C = 0.49` quoted
for the spherical gNFW² profile. Script:
`analysis/code/check_gnfw2_quadrature.py` (seconds, few hundred MB, safe to
run alongside the gce-fisher injections). Everything below is measured, not
argued.

1. **The number is under-quadratured.** The correct value at the smoothing it
   is quoted at is **0.51**, not 0.49.
2. **The number is not measured with the beam the figure paragraph
   describes.** Under the figure's own 3° top-hat the same profile gives
   **0.45**. Whichever value is quoted, the text must name the kernel.

Both move in the conservative direction — the spherical cusp is *more*
concentrated than the manuscript currently claims, so every comparative
statement strengthens. Nothing here weakens the Comment.

## What C = 0.49 is, exactly

Emitted by `stage4_stats.py:94-95` as `conc_3_10`:

> C ≡ F(<3°)/F(<10°), where F(<θ) is the **sum of map pixel values** over
> pixels with `hypot(l, b) < θ` on the 40°×40°, 0.1° grid, after **Gaussian
> smoothing at FWHM 2°**, measured on the **216-view Stage 3 suite**
> (`stage3_maps/`, 36 azimuths × 6 halos) — *not* on the figure's maps.

The gNFW² entry is the analytic template from `stage4_templates.py:129`,
`dm_templates.j_map` at `A2020_DEFAULTS`:

| | |
|---|---|
| functional | **J = ∫ρ²ds** — a genuine density-*squared* LOS integral, `trapezoid(rho**2 * jac, ...)`, `dm_templates.py:277` |
| γ | 1.2 |
| ρ⊙ | 0.28 GeV cm⁻³ |
| r_s | 26 kpc |
| R⊙ | 8.18 kpc |
| geometry | spherical, q_b = q_c = 1 |
| s_max | 100 kpc |
| LOS nodes | sinh-substituted, n_s = 800 |

**The LOS integral is not the problem.** It is converged to 3×10⁻⁶ between
n_s = 400 and 800 and to 9×10⁻⁷ between 800 and 3200. The `s = s₀ + p·sinh u`
substitution does exactly what it was written to do.

## Finding 1 — the sky-plane quadrature

`stage4_templates.py` evaluates `j_map` at pixel **centres**. That is a
midpoint rule on a γ = 1.2 cusp, which is the failure mode
`j_factor_roi`'s own docstring warns about for Cartesian grids — and it is
load-bearing here because the cusp is where the flux is:

- **14.4%** of the unsmoothed flux inside 3° lies inside 0.1°;
- the central pixel is evaluated **1.706× too low**;
- total F(<10°) is 3.3% low.

Numerator and denominator are both suppressed by nearly the same *absolute*
amount, so the ratio is biased low. Supersampling each cell (9×9 globally,
61×61 inside 2°, 401×401 inside 0.35°) gives the correct cell average:

| smoothing | pixel-centre (as published) | pixel-averaged (correct) |
|---|---|---|
| Gaussian FWHM 2° | **0.4910** | **0.5079** |
| Gaussian FWHM 3° | 0.4535 | 0.4696 |
| 3° flat top-hat | 0.4346 | 0.4530 |

The pixel-centre column reproduces the published 0.491 / 0.4535 exactly,
which is the provenance check: same statistic, same grid, same template, only
the cell quadrature changed.

**Only the gNFW² column moves.** F98 / Coleman20 / NB come from
flux-conserving FITS regrids — already cell averages — and are not cuspy;
Cao13 is gcepy's energy-summed template; the corrected J and old-star maps are
particle histograms. No other entry in the concentration table is affected.

### Downstream

| statement | was | becomes |
|---|---|---|
| gNFW² C, FWHM 2° | 0.49 | **0.51** |
| gNFW² C, FWHM 3° | 0.45 | **0.47** |
| views with C above gNFW², FWHM 2° | 6 / 216 | **0 / 216** |
| views with C above gNFW², FWHM 3° | 0 / 216 | 0 / 216 |

The corrected-J median (0.427), its range (0.34–0.50, max 0.496), the bulge
templates (0.27 / 0.27 / 0.29) and the old stars (0.231) are all unchanged.
The six views that appeared to exceed the spherical cusp at FWHM 2° no longer
do — the "6/216 (0/216)" in
`muru_quantitative_statements_2026-08-31.md:58` becomes **0/216 at both
smoothings**, and the sentence "sits near the cuspy spherical profile" can be
stated without the exception.

## Finding 2 — three different beams, three different answers

The figure paragraph says the last two columns are "smoothed with the same
3° beam". That beam is a **flat-(l, b) top-hat disk** — `figure_render.py:178`,
`kern = hypot(...) <= 3.0`, convolved onto `jmap_muru_geom` and `stars_col` at
lines 227–228 — chosen to match their aperture. It is a genuinely different
kernel from a Gaussian FWHM 3°: flatter core, hard edge. It is also applied to
a **different map product**: `figure_maps.npz`, six halos at the published
viewpoint, versus the 216-view `stage3_maps` suite the statistics run on.

So an unqualified "C = 0.49" sitting near that paragraph invites the reader to
attribute it to the 3° beam. Under that beam the same profile gives **0.435**
as published, **0.453** corrected — and 0.435 lands on top of the corrected-J
median of 0.427. That is not a real degeneracy, since the two are not measured
with the same kernel, but it is precisely the collision an unnamed kernel
invites: a referee who recomputes on the figure's footing gets ~0.43 for both.

**Recommendation: name the kernel every time C appears**, and say the
concentrations come from the 216-view suite rather than from the figure.

## Suggested wording

⚠ Not applied — Overleaf is authoritative and collaborators are editing it.
The local `comment_muru.tex` is stale and contains no C sentence at all, so
this is drafted from scratch rather than as a diff; match it to whatever the
current Overleaf sentence says.

Main text, replacing wherever 0.49 is quoted:

```latex
Smoothing every map to a common $2^{\circ}$ Gaussian beam and forming the
concentration $C\equiv F(<3^{\circ})/F(<10^{\circ})$, the corrected
annihilation maps give $C=0.43$ over the 216 viewing geometries, against
$0.27$, $0.27$ and $0.29$ for the Freudenreich, Cao and Coleman bulge
templates and $0.23$ for the same simulations' old stellar populations.
Every one of the 216 maps is more concentrated than every bulge template,
and none reaches the $C=0.51$ of a spherical gNFW$^{2}$ cusp with
$\gamma=1.2$.
```

If the sentence must instead sit on the figure's footing, the kernel changes
and so does every number in it — the whole table would need recomputing under
the top-hat, which has not been done for anything but the gNFW² template.
**Do not mix the two.**

A footnote or SM sentence worth carrying:

```latex
$F(<\theta)$ sums surface brightness over pixels within $\theta$ of the
Galactic center in flat $(l,b)$; the analytic gNFW$^{2}$ template is
evaluated as a cell average rather than at pixel centers, which matters at
the $3\%$ level for a $\gamma=1.2$ cusp on the $0.1^{\circ}$ grid.
```

## Finding 3 (partial) — the SM correlation table also moves

Re-measured on **G1.1 only** (36 views, so the shift is characterized but the
216-view medians are not):

| | median, centre → averaged | max shift |
|---|---|---|
| corr(J, gNFW²), FWHM 2° | 0.8805 → 0.8647 | 0.017 |
| corr(stars, gNFW²), FWHM 2° | 0.7295 → 0.7092 | 0.020 |
| corr(J, gNFW²), FWHM 3° | 0.9337 → 0.9254 | 0.009 |
| corr(stars, gNFW²), FWHM 3° | 0.8012 → 0.7882 | 0.013 |

(G1.1's absolute medians sit below the published all-halo values; only the
*shift* is the measurement here.)

Applied to `muru_sm_correlation_table.tex`, the gNFW² row plausibly goes
0.92 → **0.90** ($J$) and 0.72 → **0.71** (stars). Two consequences for the
SM text, which is written around that row:

- "the Coleman template ($0.93$) lying marginally above gNFW$^{2}$ ($0.92$)"
  — the gap **widens**, so the argument strengthens.
- 🚨 "At $3^{\circ}$ smoothing the ordering reverses again (gNFW$^{2}$
  $0.96$, Coleman $0.95$)" — a $-0.009$ shift takes gNFW² to $\approx 0.95$
  and the reversal may become a **tie**. The SM paragraph leans on that
  reversal to argue the correlation is not measuring a stable property. It
  still is not, but the sentence may need rewording.

**This needs the full 216-view recomputation before the SM is finalized** —
rerun `stage4_stats.py` with the cell-averaged template. The table header
already says the values are emitted programmatically and must not be
hand-edited; that instruction stands.

## What is not affected

- The corrected J maps, old-star maps, and their statistic — all particle
  histograms, no analytic cusp to under-resolve.
- Every bulge template.
- The axis ratios and c4: gNFW² is spherical, so q30 = 1.000 and c4 = 0.000
  by construction, cell averaging or not.
- The continuum test, the frame residuals, the anchor validation.
- **Fig. 1 and Fig. 2**: neither renders the analytic template.
