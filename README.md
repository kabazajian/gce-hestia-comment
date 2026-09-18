# Analysis code for the Comment on Muru et al., PRL 135, 161005 (2025)

Code, derived data products, and the figure for

> K. N. Abazajian, J. Kumar, and O. Macias, *Comment on "Fermi-LAT Galactic
> Center Excess Morphology of Dark Matter in Simulations of the Milky Way
> Galaxy"* (2026).

Every number quoted in the Comment can be recomputed from what is in
`derived/` without a particle file. The particle data are not included and
cannot be (see **Data**), but everything needed to re-run from particles is
here once you have them.

## Layout

```
code/       the analysis, in pipeline order (table below)
docs/       the plan, the validation against Muru et al.'s own code, a record per stage
derived/    the small products every number in the Comment rests on
figures/    the figure as the code produces it, plus the SM correlation table
```

## The point, in three lines

Muru et al. square a projected mass; annihilation integrates a squared
density. The estimators differ by where the square goes and by a factor s²:

```
Muru et al. "DM":          sum_i m_i               ~  int rho     s^2 ds
Muru et al. "quadratic":  (sum_i m_i)^2            ~ [int rho     s^2 ds]^2
annihilation (this work):  sum_i m_i rho_i / s_i^2  ~  int rho^2       ds
old stars   (this work):   sum_i m_i       / s_i^2  ~  int rho_*       ds
```

`code/particle_projection.py` is the estimator library and
`code/test_particle_projection.py` its tests (`pytest`; numpy and scipy only).

## Pipeline

| step | script | what it does |
|---|---|---|
| 0 | `anchor/anchor_muru_g32.jl`, `anchor/anchor_compare.py` | runs Muru et al.'s own Julia code on one halo and checks our implementation of their statistic against it, sight line by sight line: **5×10⁻⁷** over all 160,801 sight lines of their grid |
| 1 | `stage1_frames.py` | per-halo center, disk plane, bar angle and A₂ |
| 2 | `stage2_density.py` | k-nearest-neighbour density at each DM particle (k = 32; 16 and 64 checked) |
| 3 | `stage3_views.py` | 216 views (36 azimuths × 6 halos), four maps each |
| 4 | `stage4_templates.py`, `stage4_stats.py`, `stage4b_continuum.py` | comparison templates, morphology statistics, hold-out test |
| 5 | `stage5_systematics.py`, `stage5_summary.py` | systematics budget |
| — | `ahf_frames.py`, `ahf_azimuth_scan.py`, `ahf_radius_scan.py` | viewing frames from the AHF profiles; see `docs/viewing_frames.md` |
| — | `check_gnfw2_quadrature.py`, `gnfw_jmap.py` | the gNFW² comparison template and its pixel-quadrature check (`C = 0.51`) |
| — | `conc_per_halo.py` | per-halo concentrations and the unresolved-centre test (SM Table I) |
| — | `figure_maps.py`, `figure_render.py` | the figure |

Stages 1–3, 5 and `figure_maps.py` need the particle files. Everything else
runs from `derived/`.

### Reproducing the figure without particles

```
cd code
python figure_render.py --all --tag bestaz     # reads ../derived/figure_maps.npz
```

`figure_render.py` expects the maps beside the data; point `NFS` at
`../derived` or symlink. It prints the contour levels and axis ratios
against the values printed in Muru et al.'s Figs. 1–2.

## Validation, in short

* Our implementation of Muru et al.'s projection statistic reproduces their
  published Julia code, run on the same particles in a common geometry, to
  5×10⁻⁷ (`docs/muru_anchor_2026-08-30.md`). The residual is their Float32
  accumulator.
* The corrected estimator was validated against the analytic ∫ρ²ds of a
  sampled NFW profile before the simulation data arrived.
* With the viewing frames from the AHF profiles Muru provided, the published
  axis ratios are reproduced to 0.02 (DM) and 0.05 (stars) — details and two
  traps in `docs/viewing_frames.md`.
* No conclusion moves under kNN k, line-of-sight truncation, observer
  distance, particle half-sampling, centering, smoothing width, or excision
  of the unresolved centre (`docs/muru_stage5_systematics_2026-08-30.md`,
  `docs/conc_table_suggestion_2026-09-16.md`).

## Data — not included, and why

The six HESTIA particle files (4.6 GB) and the AHF profile files are
**M. M. Muru's to distribute**. The PRL data-availability statement says the
simulation data are "available from the authors upon reasonable request",
which is how they reached us. Please request them from him; this repository
is not a redistribution channel. Two large intermediates are also omitted
because they regenerate in minutes: `stage3_maps/` (449 MB) and `rho_k32/`
(172 MB).

What is here is everything the Comment's numbers rest on:

| file | contents |
|---|---|
| `derived/figure_maps.npz` + `_meta.json` | the maps behind the figure, all six halos, at the frames used |
| `derived/stage4_stats.json` | every morphology statistic, all 216 views, both smoothings |
| `derived/stage5_systematics.json` | the systematics variants |
| `derived/stage4_templates.npz` | the bulge and gNFW² comparison templates (gNFW² cell-averaged) |
| `derived/conc_per_halo.json` | per-halo concentrations and the core-excision test |
| `derived/muru_frames.json`, `ahf_frames*.json` | per-halo centers, disk normals, bar axes; the AHF-derived viewing frames |
| `derived/ahf_*_scan.json` | the radius and azimuth scans of the viewing frame |
| `derived/figure_ratios.json` | the contour levels and axis ratios printed on the figure |

One filter is applied: statistics against `combo_*` templates — bulge +
diffuse combinations weighted by fitted amplitudes from a separate,
unpublished analysis — are removed. The Comment uses only the individual
bulge templates and gNFW², so nothing it says depends on them.

## Dependencies

`numpy`, `scipy`, `matplotlib` (with a LaTeX installation for `text.usetex`),
and `pyarrow` to read the particle files. The gNFW² line-of-sight integrator
is vendored in `code/gnfw_jmap.py`. **Template construction**
(`stage4_templates.py`) additionally needs the public bulge-template FITS
files, the `gcepy` package, and the unpublished fits; it is included for the
record, and the constructed templates are shipped in `derived/` so nothing
downstream requires it. Stage 0 needs Julia 1.11 and Muru et al.'s
repository.

## Provenance and attribution

* Muru et al.'s analysis code is theirs, at
  <https://gitlab.aip.de/muru/gce_in_hestia>; it is referenced, not
  bundled. `code/anchor/HestiaUtils/` is a stub we wrote so their sources
  load without their private data loader — every function in it raises.
* Bulge templates, redistributed here only as regridded comparison maps:
  the Freudenreich (1998) and Coleman et al. (2020) maps from C. Gordon's
  public repository, <https://github.com/chrisgordon1/galactic_bulge_templates>;
  the nuclear bulge from the template release accompanying Pohl, Macias,
  Coleman and Gordon (2022); and the Cao et al. (2013) template from `gcepy`
  (S. D. McDermott et al.), <https://github.com/samueldmcdermott/gcepy>.
* Portions of this code were written with Claude Code (Anthropic) under the
  authors' direction, as stated in the Comment's acknowledgments.
