# Muru et al. (arXiv:2508.06314): did they compute the annihilation projection correctly?

**The point is due to Jason Kumar**, raised in a journal club and relayed by Kev,
2026-07-30. The annihilation signal is ∫ρ² dℓ along the line of sight. Did they
instead square the projected ρ? And is "indistinguishable from stars" credible
when DM enters as ρ² and stars as ρ?

⚠ **Attribution matters here and an earlier draft of this file got it wrong**,
crediting the question to Kev. It is **Jason Kumar's** point. If this becomes a
Comment (see `validation/muru_comment_prep.md`), that ordering is not a
courtesy — it is authorship.

**Answer: no, they did not compute it correctly, and the error runs in exactly
the direction that manufactures their conclusion.** Confirmed twice over — from
their Appendix A text *and* from their own published analysis code. Two independent problems,
plus a resolution limit that would bite even if both were fixed.

---

## 0. Venue

arXiv metadata records: *"Accepted for publication in **Physical Review
Letters**, 9 pages, 3 figures, 1 table"*. There is **no `journal_ref` or DOI**
on arXiv, and the PDF in `docs/` is **v1 (2025-08-08)**.

## 1. What they actually did — their Appendix A, verbatim

> **4.** We filter the particles for the suitable particle type … and remove
> particles further than 15 kpc from the observer. The particles are then
> projected on a unit sphere using galactic coordinates.
>
> **5.** We create a regular grid of sight lines … Each sight line (l, b) is
> assigned a density value by **summing up the masses of all the particles
> within 3° of the sight line**.
>
> **6.** For quadratic density, **the previously calculated DM density values
> are squared**.
>
> **7.** All the density values are normalized using the maximum density value.

The "previously calculated DM density values" of step 6 are the step-5
quantities, which are already aggregated along the line of sight. So they form

  **their map ∝ [ Σ masses in a 3° cone ]²  ∝  [ ∫ ρ(s) s² ds ]²**

against the physically required

  **J(ℓ,b) = ∫ ρ²(s) ds**

### Problem A — squaring the projection instead of projecting the square

`(∫ρ ds)² ≠ ∫ρ² ds`. This is Jason Kumar's point, confirmed from their own text
and, decisively, from their own code (§1a).

### Problem B — mass-in-a-cone is not even a column density

Summing particle *masses* within an angular aperture gives ∫ρ dV = ΔΩ ∫ρ s² ds
— an **s²-weighted** integral, because the cone's volume grows with distance.
The physical J-factor carries **no** s² factor (flux dilution 1/s² exactly
cancels the volume s²). So their pre-squaring quantity over-weights distant
material; combined with the 15 kpc truncation this is an ad hoc radial window,
not a line-of-sight integral.

## 1a. 🚨 CONFIRMED IN THEIR OWN PUBLISHED CODE

Their ref. [54], `https://gitlab.aip.de/muru/gce_in_hestia`, cloned 2026-07-30.
The text could conceivably have been loose wording; the code cannot be.

**The aggregation function is named and documented as a mass sum**
(`src/BulgeDensityProjection.jl`):

```julia
"""Calculate the line-of-sight mass for each `l`, `b` pair in `grid` ...
Sums together all the masses within `distance` from the line-of-sight."""
function calc_density_on_grid(particles, tree, grid, distance)
    masses = particles.Masses
    los_mass = zeros(Float32, size(grid)...)
    for i in eachindex(grid)
        los_mass[i] = inrange(tree, grid[i], distance) |> (idxs -> masses[idxs]) |> copy |> sum
    end
    los_mass
end
```

The variable is literally `los_mass`, and the KDTree is built on `(l, b)` only —
so `inrange` is an **angular** neighbour search and the result is the total mass
in a 3° aperture, integrated to the 15 kpc cut. Confirms Problem B.

**The squaring is applied to that grid, on the exact path that makes the paper's
figures:**

```julia
# src/plot_methods_densityprojections.jl:349
dens = Dict(:dm => dens_dm[i], :stars => dens_stars[i], :dm2 => dens_dm[i].^2)

# scripts/plot_6_galaxy_projections.jl:81  (their Figs. 1-2)
plot_6_density_projections(map(d -> d.^2, dens6_dmonly)[dmonly_order], ...)

# src/calc_methods_densityprojections.jl:116
dens = calc_density_on_grid(particles, lbs_grid, max_dist_from_los) .^ density_exponent
```

`:dm2` — the "quadratic DM" panel — is `dens_dm .^ 2`, the element-wise square of
the already-aggregated line-of-sight mass map. Confirms Problem A. There is no
path in the repository that squares the density *before* aggregating.

**What this does NOT overturn, stated so the Comment is fair.** Their NFW
reference is computed the same way (`scripts/plot_nfw.jl:33`,
`plot_density_projection(dens_mock.^2; ...)`), so their *internal* comparison of
simulated halos against a spherical NFW uses one consistent statistic and the
relative asphericity claim is unaffected by this error. What fails is the step
that matters: comparing that quantity against a **linear** stellar projection and
concluding the γ-ray morphologies are indistinguishable.

## 2. Why this is not a small error: squaring is a *relabeling*

**Squaring is a monotonic pointwise transform, so it leaves every isophote
shape unchanged and only renames the levels.** The contour of f² at level c² is
the identical curve to the contour of f at level c.

Verified numerically on a triaxial gNFW (γ = 1, q_y = 0.85, q_z = 0.70,
observer 8 kpc in the plane):

```
axis ratio of [∫ρ ds]²  at 50%  of its max   = 0.6667
axis ratio of  ∫ρ ds    at 70.7% (=√0.5)     = 0.6667      identical
```

So their "quadratic DM" panels carry **no morphological information beyond the
linear DM panels**. Whatever shape they report for the annihilation map is the
shape of the *column density* map, sampled at a different contour level.

That is the mechanism behind the result Jason Kumar found incredible. Correctly, DM
enters as a ρ²-weighted projection and stars as a ρ-weighted one — genuinely
different functionals. Under their procedure DM becomes a monotonic function of
a ρ-weighted projection, i.e. **the same shape family as the stars**. The
"indistinguishable" conclusion is substantially an artifact of the method.

## 3. How different is the correct calculation? Dramatically.

Same model, half-width along ℓ at each contour level:

| level | ∫ρ ds (column) | their (cone-mass)² | **correct ∫ρ² ds** |
|---|---|---|---|
| 70% | 0.15° | 0.05° | ≲0.01° |
| 50% | 0.90° | 0.15° | ≲0.01° |
| 30% | 4.60° | 0.65° | **≲0.05°** |
| 10% | 20° | 4.55° | 0.15° |

The ρ² weighting concentrates the signal by **one to two orders of magnitude in
angular scale**. A correct annihilation map is a compact, strongly peaked
structure; a column-density map is an extended one. They are not the same object
and a template fit distinguishes them easily.

(The exact numbers are model-dependent — a cuspy γ = 1 profile maximizes the
contrast — but the direction and rough size are robust to profile choice.)

## 4. A separate problem that survives fixing both: resolution

| quantity | value |
|---|---|
| HESTIA softening length | 0.22 kpc |
| subtended at 8 kpc | **1.58°** |
| their smoothing aperture | **3.0°** |
| half-width of correct ∫ρ² ds at 30% of max | ≲0.05° |

**The region that carries essentially all of the annihilation signal lies inside
the softening length and far inside their 3° aperture.** So even with the
integral done correctly, these simulations could not resolve the annihilation
morphology where it matters. This is a limitation of resolution, not of care,
and it is worth stating separately from the two errors.

## 5. How much would a *correct* map match our bulge templates? (qualitative)

Qualitatively, **much less well than they claim**, for two reasons that act
together:

1. **Concentration.** The Coleman20 / F98 boxy bulge templates are *extended*
   structures spanning several to ~10°, tracing the stellar bar, with genuine
   boxy (non-elliptical) isophotes. A correct ρ²-weighted annihilation map is
   compact and strongly centrally peaked. Even if the *shapes* of the isophotes
   were similar, the radial profiles are not, and a bin-by-bin template fit
   keys on exactly that.
2. **Which part of the halo sets the shape.** ρ² weights the innermost, densest
   material, so the annihilation morphology is set by the *inner* halo. Their
   own result is that flattening *increases outward*, so weighting inward should
   make the correct map **rounder** than what they display — the opposite of
   their claim that squaring "amplifies" the anisotropy.

**The honest caveat:** the inner halo need not be round. A bar and disc can
flatten and triaxialize the central DM, and that is a real effect their paper is
right to raise even if their demonstration does not establish it. So the correct
statement is not "DM is spherical" but "**their analysis does not measure the
annihilation morphology**, and the physically-correct version would be more
centrally concentrated and plausibly rounder than what they show."

## 6. What a quantitative answer needs

Their density profiles, as Kev noted. Two routes:

- **Their public repository**, cited as their ref. [54]:
  `https://gitlab.aip.de/muru/gce_in_hestia`. This is the clean route — it
  should permit recomputing ∫ρ² ds correctly from the same halos and comparing
  directly against our bulge templates.
- Failing that, digitizing the axis ratios printed in their Figs. 1–2 — but note
  those are the *relabeled column-density* shapes (§2), so they cannot be
  reinterpreted as annihilation morphologies without the underlying profiles.

⚠ Per rule 1a, nothing quantitative about their halos enters `validation/`
without one of those sources. Everything in this document is either quoted from
their PDF or computed here from a stated toy model.

## 7. Consequences for our own analysis

The Task 9 **boxy-DM injection test** added earlier today is still worth doing —
but its motivation changes. It is no longer "Muru et al. show DM is boxy, check
we are not fooled"; it is "**morphological degeneracy between a DM template and
the stellar bulge is a real systematic direction, whoever raises it**". The test
should span a range of DM morphologies from round to flattened-and-boxy, and
report where the bulge starts absorbing signal — which is a statement about our
limit's robustness, independent of whether Muru et al.'s specific claim holds.
