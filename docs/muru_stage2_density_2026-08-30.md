# Muru redo Stage 2 — kNN density: k = 32 adopted, caches built

**2026-08-30.** Plan: `validation/muru_redo_plan_2026-08-30.md` §2. Script:
`scripts/muru_redo/stage2_density.py`, run on the Studio against the NFS
mount of nucosmo's data (`/private/nfs/Data2/hestia_muru`).

## k-convergence on G3.2 (the anchored halo), k ∈ {16, 32, 64}

**Per-particle ρ_k/ρ_32** — median within 0.3% of unity in *every* radial bin
from 0.22 kpc (the softening) to 30 kpc: the estimator has no k-dependent
bias anywhere we will use it. The 16–84% spread (±12% at k=64, ±21% at k=16)
is per-particle scatter, which the J integral averages down by design.

**The quantity actually consumed** — correct J-maps from one solar-circle
view (R⊙ = 8.18 kpc along the bar axis, ±20° grid, σ = 1° Gaussian smoothing),
relative to k = 32 over pixels above 1% of max:

| k | median | max |
|---|---|---|
| 16 | 3.3% | 12.2% |
| 64 | **1.6%** | **4.3%** |

**Adopted: k = 32**, carrying **~2% (median) – 4% (max) at σ = 1°** as the
density-estimator systematic. Stage 5 re-measures it at the actual Stage 4
smoothing scales (2–3°), where it can only shrink.

## Units / normalization validation on real data

Median kNN ρ against the direct shell M/V on G3.2: ratio 1.03 (0.5–1 kpc)
rising smoothly to 1.19 (32–64 kpc). No constant factor — an h slip or a
4π/3 error would sit at 1.5–2× — and the absolute scale is sane (ρ at
8–16 kpc ≈ 4.8×10⁶ M☉/kpc³ ≈ 0.18 GeV cm⁻³ for this M ≈ 0.9×10¹² M☉ halo).
The mild outward-growing excess is the expected mass-weighted sampling of an
aspherical, substructured halo (the median particle sits in denser-than-
shell-average regions), not an estimator defect — and it is common to all k,
so it cancels in the convergence numbers above.

## Production caches

`nucosmo:/Volumes/Data2/hestia_muru/rho_k32/` (= NFS
`/private/nfs/Data2/hestia_muru/rho_k32/`): one float32 `.npy` per halo, DM
rows in file order, physical M☉/kpc³, ~45M values total, with
`PROVENANCE.json` (script, k, date, per-halo n and ρ stats). Regenerable in
~4.5 min; not in git.

`knn_density` gained chunked queries (`block=2_000_000`) — the all-at-once
query on a 10M-particle halo allocated ~10 GB; identical results, tests pass
(5/5), peak memory now ~2 GB per halo.
