# Stage 3a anchor: our `project_muru` vs Muru et al.'s own pipeline — PASS

**2026-08-30.** Plan step: `validation/muru_redo_plan_2026-08-30.md` §3a. The
question the anchor answers: does `gcereduce.particle_projection.project_muru`
reproduce the *actual published pipeline* of arXiv:2508.06314 — frames, units,
handedness, aperture — on real HESTIA particles, before any production map is
made?

**Answer: yes, to 5×10⁻⁷, on all 160,801 sightlines of their exact figure
grid — with one convention discovered and adopted along the way (flat-(l,b)
aperture metric), and one bookkeeping error caught by the discrimination
layer (an ℓ-mirror wrongly applied at first, failing at 26%).**

## Setup

- **Their code, verbatim.** `gitlab.aip.de/muru/gce_in_hestia` re-cloned to
  `~/work/gce_in_hestia`, HEAD `b50114d9` (2026-01-23) — **no commits since
  the 2026-07-30 analysis clone**, so the state analyzed in
  `validation/muru_projection_analysis.md` is the current state (Comment
  item 4.2 partly discharged). ⚠ Item 4.2 note: the HEAD commit *"feat: add
  possibility to square density"* is dated 2026-01-23, i.e. **after** arXiv
  v1 (2025-08-08) — the public repo reflects a post-submission state; worth
  one line in the Comment file, not over-interpretation.
- **Private dependency stubbed, not faked.** Their `HestiaUtils` is a local
  package on their machine (`path = /z/moortis/.julia/dev/HestiaUtils`; rule
  3). `scripts/muru_redo/anchor/HestiaUtils/` is a stub with the same UUID
  whose every loader **errors if called** — their two src modules
  (`GalacticCoordinates.jl`, `BulgeDensityProjection.jl`) load unmodified,
  and geometry is supplied explicitly instead of through their AHF readers.
- **Julia 1.11.4** (their Manifest's version) via juliaup on the Studio.
  nucosmo's macOS 10.13 cannot run Julia 1.11; the G3.2 file (smallest halo,
  474 MB) was copied to `~/work/hestia_muru/` on the Studio. ~10 s runtime.
- **Appendix A parameters, their hard-coded values**: observer 8·0.677 kpc/h
  from the center (α = 0 from the supplied disk vector), particle cut
  d ≤ 15·0.677 kpc/h, aperture 3.0°, grid −20:0.1:20 (401×401). Frame
  vectors deliberately **non-axis-aligned** so a handedness error cannot
  hide. Working units are kpc/h throughout (their `coordsinMpc` ×10³ with
  centers pre-scaled — the center passed to their code must be kpc/h).

## The three layers (`scripts/muru_redo/anchor/anchor_compare.py`)

| layer | measured | gate |
|---|---|---|
| per-particle (l, b, d) transform, 200k particles | max abs diff 8.5×10⁻¹⁴ / 3.3×10⁻¹³ / 3.6×10⁻¹⁵ | < 10⁻⁴ |
| map, `project_muru(metric="flat_lb", exponent=1)`, 160,801 cells | max rel 5.0×10⁻⁷, median 1.8×10⁻⁷ | < 10⁻⁵ |
| discrimination: ℓ-mirrored | max rel **0.28** | must exceed 0.1 |
| discrimination: `metric="flat_lcosb"` | max rel **0.10** | must exceed 10⁻³ |

The 5×10⁻⁷ residual is their Float32 accumulator (`los_mass = zeros(Float32,…)`)
against our float64 — exactly the expected scale, no other source.

## Two findings worth keeping

1. **Their 3° aperture is flat in raw (l, b) degrees — no cos b.** Their
   KDTree is built on (l, b) directly (`calc_KDTree`), so the search disk is
   not a constant angular aperture: in true angle it widens as 1/cos b, ~6%
   elongated at |b| = 20°, worth up to **10% in summed mass** on the G3.2
   map. `project_muru` gained `metric="flat_lb"` to reproduce this exactly;
   the correct-J estimator keeps the true small-angle metric. Scope
   discipline (Comment prep §5): both their DM and stellar maps go through
   the same aperture, so this is a *description accuracy* detail, not a new
   objection — record it, do not add it to the Comment's claim.
2. **The frame mirror is real but belongs to frame construction, not to the
   transform.** Their +ℓ axis is e_gc × e_n where our library's is
   e_n × e_gc — a first comparison applied that mirror on top of pre-rotated
   coordinates (where both codes share axes by construction) and failed at
   26%. The discrimination layer is what caught it. Production Stage 3 maps
   built in our own frame and compared against anything of theirs must
   handle this ℓ-mirror — the project's eighth handedness trap, recorded
   before it cost anything.

## What this validates, and what it does not

Validated: `project_muru(metric="flat_lb")` **is** their pipeline (transform,
units, aperture, statistic) on real particles. The Stage 3 reproduction leg
can be run entirely in Python.

NOT validated by this anchor: the correct-J estimator (that validation is
analytic, `tests/test_particle_projection.py`), their *published figure
values* (needs their actual AHF axes — requested from Muru), and anything
about the stellar maps (their star filter `0 < a_form < 0.79`, i.e. old
stars only, is noted for the Stage 3 stellar leg).

Tests: `tests/test_particle_projection.py` — 5 pass (new
`test_flat_lb_metric_matches_their_aperture_not_ours` keeps the two metric
branches from being conflated in refactors).
