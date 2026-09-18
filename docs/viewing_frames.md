# Viewing frames: how the columns of Fig. 1 that reproduce Muru et al. are oriented

Appendix A of Muru et al. (PRL 135, 161005) places the observer 8 kpc from
the center along the major axis of the inertia tensor of particles within
10 h⁻¹ kpc, with the minor axis as the disk normal. Their code takes these
axes from the AHF halo-finder profile file (`get_galvectors_from_profile`,
in a package not included in their public release). M. M. Muru provided the
profile files; what follows is what they contain and how we use them.

## The disk normal is well determined; the in-plane major axis is not

At the profile bin nearest 10 kpc/h the six halos have b/a = 0.977–0.994.
The major and intermediate axes are therefore within 1–2% of each other, and
the direction called "major" in the plane is poorly conditioned: it rotates
by 25°–82° between adjacent radial bins of the same profile, while the minor
axis (c/a ≈ 0.7) moves by less than 1.5°. The minor axis agrees with the
stellar disk normal computed independently from the particles to 0.1°–2°.

A scan of the observer azimuth through the full plane
(`code/ahf_azimuth_scan.py`) shows the stellar projection is a strong,
180°-periodic function of azimuth — the bar — while the DM projection varies
much less. The published stellar and DM panels together therefore pin the
azimuth to within about 10°.

## What the figure uses

The disk normal is the AHF minor axis. The in-plane azimuth is chosen per
galaxy, within the degeneracy above, to match the published corner tables
(`derived/ahf_frames_bestaz.json`; azimuths lie 10°–80° from the major axis
of the 10 kpc/h bin). The published axis ratios are then reproduced to 0.02
(DM) and 0.05 (stars), and the published contour levels to within six
percentage points. With the observer placed exactly on the major axis of
that bin instead (`derived/ahf_frames.json`), the residuals are 0.03–0.14
(DM) and 0.05–0.22 (stars) — the difference is the degeneracy, not the
pipeline: our implementation of their projection agrees with their public
code to 5×10⁻⁷ on all 160,801 sight lines of their grid given a common
geometry (`docs/muru_anchor_2026-08-30.md`).

## Two details anyone re-running this should know

**Halo labels.** The particle files as distributed are named G11, G12, …,
and these names are transposed within each pair relative to the paper:
particle file G11 is the paper's G1.2 and G12 is G1.1, and likewise for the
G2 and G3 pairs. This was established two ways — the AHF profile files match
Table I of the paper by both M200c and R200c under their own names, and each
profile matches the *other* particle file of its pair on the cumulative mass
profile to 0.2% — and it is confirmed by the fact that the published panels
are reproducible only under the swap. Every key in this repository's data
products is the particle-file label; `figure_render.PAPER_LABEL` maps to the
paper's names at presentation time. Do not rename the data keys.

**The b = 0 level rule.** The paper matches DM contour levels to the stellar
contours "at the same l coordinates" where they cross b = 0. A stellar
contour crosses b = 0 at two longitudes, which differ for a barred map, and
the paper does not say which is used; the public code selects, per level,
whichever b = 0 pixel has the stellar value nearest that level, on either
side of the center. We use the positive-ℓ crossing throughout
(`figure_render.matched_levels`). The two choices differ by up to five
percentage points in the lowest DM level, and applying the public code's
rule to our maps changes no axis ratio by more than 0.04.
