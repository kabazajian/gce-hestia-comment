# Muru et al. (arXiv:2508.06314) corrected redo on the HESTIA particles — plan

**Approved by Kev 2026-08-30** (decisions §7). This is Comment-prep item 4.3
(`validation/muru_comment_prep.md` — "do the corrected calculation ourselves,
before submitting") extended into a distribution over viewing geometries, per
Kev's request of the same date.

**Goal.** Compute the correct annihilation morphology J(ℓ,b) = ∫ρ²ds for the
six HESTIA halos over many observer positions on the solar circle; build the
empirical distribution of simulated GCE morphologies; test whether that
distribution is a describable continuum (hold-out validation); and locate the
data-preferred bulge model(s) within it. Muru et al. squared the projected
(line-of-sight–aggregated) density — `[∫ρ s² ds]²` — where the annihilation
signal is `∫ρ² ds`; the full diagnosis is
`validation/muru_projection_analysis.md`.

**Data.** `nucosmo:/Volumes/Data2/hestia_muru/` — six
`HESTIA_*_minimal_particles.arrow` files (G1.1/G1.2, G2.1/G2.2, G3.1/G3.2 =
3 Local Group realizations × 2 galaxies), from M. M. Muru via the AIP share,
2026-08-30, each verified byte-exact. Columns: ParticleIDs, ptype
(dm/bh/gas/stars), Coordinates1-3 (absolute box, **Mpc/h, h = 0.677**),
Velocities1-3, Masses, GFM_StellarFormationTime (scale factor; ≤ 0 ⇒ stellar
wind, not a star). AHF halo centers are in the README beside the data.

**Compute.** Entirely on nucosmo (`micromamba run -n hestia`, data local).
Zero load on the Studio's Task 9 queue. Scripts in `scripts/muru_redo/`,
results to `results/`, run records under `runs/` as usual.

**Estimator.** `src/gcereduce/particle_projection.py`, validated in advance
against an analytic NFW (`tests/test_particle_projection.py`, 4 pass):
`project_annihilation` implements J(Ω)·ΔΩ = Σ_{i∈beam} m_i ρ_i / s_i² with
ρ_i from kNN; `project_muru` reproduces their statistic on identical
particles. One extension needed: a stellar flux map, weight m_i / s_i²
(same `_beam_sum` core), with its own test.

---

## Stage 0 — provenance and unit checks

- Convert box coordinates to halo-centered physical kpc:
  `(x − c_AHF) × 1000 / 0.677` (z = 0, comoving = physical).
- Checkable before any science: G_n.1–G_n.2 pair separations must be Local
  Group scale (~0.5–1 Mpc; from the centers alone: 0.87 / 0.68 / 0.85 Mpc);
  enclosed DM mass within 250 kpc ~ 10¹² M☉ (TNG mass unit: 10¹⁰ M☉/h —
  verify against a sanity range, not a target).
- **Email Muru** (draft: `validation/muru_data_request_2026-08-30.md`): the
  AHF profile files he offered (inertia-tensor axes), the softening length of
  these runs (we cite 0.22 kpc from the paper), and the DM particle mass.
  His replies enter `validation/reference_values.yaml` as human-supplied,
  cited to the correspondence. Also progresses Comment item 4.2.

## Stage 1 — per-halo Galactic frame

- Center: AHF, refined by shrinking spheres on DM+stars. Gate: refinement
  moves < a few softenings, else stop and investigate.
- Disk plane: net angular momentum of true stars (StellarFormationTime > 0)
  at 3–15 kpc defines ẑ. Cross-checks: stellar inertia-tensor minor axis;
  AHF axes when they arrive. Large misalignment ⇒ flag the halo.
- Bar: m = 2 Fourier mode of the face-on stellar surface density → φ_bar per
  halo, so every view records its angle to the bar (MW: ~27–30°).

## Stage 2 — density estimation

- `knn_density` (k = 32) on DM particles; k ∈ {16, 32, 64} convergence check
  on one halo before production. Wind particles excluded from all stellar
  quantities; bh/gas excluded from ρ_DM.

## Stage 3 — the map suite

Observers on the solar circle at **R⊙ = 8.18 kpc**
(`reference_values.yaml::a2020_defaults`; recorded there as our choice — A2020
never states one), every 10° of azimuth → 36 views × 6 halos = **216 views**.
Per view, on the analysis grid (40°×40°, 0.1° pixels, ±20° per
CONVENTIONS.md § Analysis geometry):

1. **Correct J-map** — `project_annihilation`.
2. **Their statistic** — `project_muru` (reproduction leg only).
3. **Stellar flux map** — Σ m_i/s_i² over wind-cleaned stars, for the
   DM-vs-stars redo with the right functional on each side.

s_max: 100 kpc for the primary leg (their 15 kpc truncation only in the
reproduction leg).

### Stage 3a — the anchor (falsifiable, before production)

With **their** geometry — observer at 8.0 kpc on the major axis, 3° aperture,
15 kpc cut — our `project_muru` must reproduce **their own published Julia
pipeline** (`gitlab.aip.de/muru/gce_in_hestia`, cloned 2026-07-30) **run on
these same particle files**. Numerical agreement gates production: frames,
units, centers, handedness all anchored against their code, with no number
read off their figures (their in-figure SVD axis ratios stay qualitative
context, rules 1a/1b). Contingency if Julia won't run on nucosmo's macOS
10.13: run their code on the Studio against one halo.

⚠ Seven coordinate-convention traps on this project so far (CONVENTIONS.md
§ THE SEVENTH). Every map comparison here includes a deliberate
flip/mirror discrimination test, not just a correlation.

## Stage 4 — morphology statistics and the bulge comparison

All maps smoothed to a common resolution ≥ the softening's subtense (1.58°
at 8 kpc; statistics reported at 2° and 3° smoothing). Per map:

- correlation/cosine similarity against the template set (below);
- SVD axis ratio per isophote level (their statistic, for continuity);
- boxiness c₄ Fourier coefficient of isophotes (boxiness ≠ ellipticity);
- concentration J(<1°)/J(<10°).

**Template set** (Kev, 2026-08-30): the best-fit bulge + diffuse combos from
the main-analysis rankings — **top 5–10 cells** of
`results/gnfw2_vs_bulge_evidence*.json` by ln H, each contributing its
best-fit bulge+NB at fitted amplitudes — plus spherical gNFW² (γ = 1.2) and
the bulge/NB pieces separately. ⚠ The A42 **bulge-λ composite** is what "λ
superimposed" ultimately means and it is *absent from the rankings so far*
(`results/gnfw2_vs_bulge_ranking_summary_2026-08-23.md` final caveat); when
it enters the rankings it becomes the primary comparison template. Until
then the top single-bulge combos stand in, stated as such.

Push each template through the identical statistic code, then report where
it falls in the empirical distribution of 216 views — and the flip side:
in how many views the corrected J-map is better matched by spherical gNFW²
than by the bulge. Expectation from `muru_projection_analysis.md` §5:
corrected maps are more concentrated and rounder; any view that isn't is
the interesting one.

### Stage 4b — continuum test via structured hold-out (Kev, 2026-08-30)

Test whether the 216 views form one describable continuum: hold out ~10%,
fit a candidate distribution to the rest, check the held-out views fit.
Three rungs of increasing strictness — a random hold-out alone would pass
trivially, since 10°-adjacent views are near-duplicates:

1. random ~10% (tests azimuthal smoothness only);
2. contiguous ~36° azimuthal arcs, one per halo;
3. **leave-one-halo-out (6-fold) and leave-one-realization-out (3-fold,
   both galaxies of a pair out together)** — the load-bearing rung: it
   distinguishes "continuum" from "six halo-specific clumps".

Space: the per-view statistic vector (~6 dims). Candidate distribution:
multivariate Gaussian after per-coordinate transforms; a 2–3 component
mixture only on visible failure; nothing more flexible — a model that
cannot fail cannot validate ("if functional and parameterizable" is doing
real work). Pass/fail per fold as a counting statement: held-out typicality
ranks vs uniform, "k of n inside the fitted 90% region".

- Continuum validates ⇒ quote the bulge templates' typicality ranks under
  the fitted distribution, with the raw counting statement alongside.
- Continuum fails ⇒ per-halo counting only; the failure is itself a finding
  (halo-history-specific morphology undercuts any general "DM looks like
  the bulge" claim).

Mirror symmetry: φ and φ+180° give mirrored maps; verify the statistics'
mirror-invariance numerically and state the effective view count honestly
(~18 per halo if invariant).

**Framing rule (Ryan's standing objection, CLAUDE.md § Writing to
collaborators):** 3 independent realizations, correlated views. Counting
statements and ranges; no σ, no CI, no p-value dressed as inference over a
population of halos that does not exist. The strict-rung effective sample
size is 3 and the write-up says so.

## Stage 5 — systematics budget

Sensitivity of every Stage 4 statistic to: smoothing scale (1.58°/2°/3°),
kNN k, s_max, center choice (AHF vs shrinking spheres), particle noise
(half-sample splits), R⊙ (8.18 vs Muru's 8.0). Claims restricted to scales
above ~2–3 softenings; the inner ~1 kpc is unresolved and we say so — that
limitation survives the corrected calculation
(`muru_projection_analysis.md` §4).

## Stage 6 — outputs and venue

`results/muru_redo_<date>.md` + figures. Feeds, in order: (a) the PRL
Comment's protective calculation (Comment prep item 4.3; Jason Kumar leads,
whoever performs the corrected calculation co-authors); (b) the paper's
§12.D paragraph (`paper/WORKING_NOTES.md`); (c) possibly separate-paper
material (Kev's earlier judgement, WORKING_NOTES ~538). Which venue gets
what is Kev's call once the numbers exist.

**✅ RULED, Kev 2026-08-30 (Stages 0–5 complete):** (a) the PRL Comment
proceeds — draft in `validation/muru_comment_draft_2026-08-30.md`, Jason
Kumar leads; (b) the paper carries one paragraph — draft in WORKING_NOTES
§12.D; (c) the separate paper is deferred and **reframed as a multi-suite
study**: hydro simulations from several groups, to get variance across the
simulations' inputs, physics and numerical treatments — a requirement our
own leave-one-halo-out failure (Stage 4b) already demonstrates within one
suite. Remaining externals: published-PRL text check, Muru's AHF axes
(exact figure frames), Jason Kumar coordination.

## Boundaries

- **Blinding untouched.** No Fermi data appears; template-vs-simulation
  only. Bulge amplitudes come from already-permitted (6a) fits; no new fit
  to real data occurs here.
- **A45/registry not implicated** — nothing here injects. If a
  HESTIA-derived template later enters Task 9's injection suite it needs
  its own registry row and PI ruling; the `gnfw/boxy` row is SET ASIDE
  (Kev 2026-08-29) and this work must not revive it by side door.
- Rule 1a: no number from their figures; the anchor is their code run on
  the same particles. Rule 1b: their printed axis ratios are context only.

## §7 Decisions recorded (Kev, 2026-08-30)

1. **λ-superimposed bulge** = best-fit bulge + diffuse combo from the main
   analysis rankings; test the top 5–10 combos. The A42 λ-composite is the
   eventual referent and is not yet ranked — tracked above.
2. **Conventions**: CONVENTIONS.md values for the primary leg (R⊙ = 8.18,
   40°×40° at 0.1°); Muru's 8.0 kpc / 3° / 15 kpc only for the
   reproduction anchor.
3. **Data request to Muru**: send now (draft file above; Kev sends).
