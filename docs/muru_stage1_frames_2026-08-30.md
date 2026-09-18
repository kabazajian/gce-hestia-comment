# Muru redo Stage 1 — per-halo Galactic frames: COMPLETE, no flags

**2026-08-30.** Plan: `validation/muru_redo_plan_2026-08-30.md` §1. Script:
`scripts/muru_redo/stage1_frames.py`, run on nucosmo; machine-readable frames
in `results/muru_frames_2026-08-30.json` (also beside the data as
`nucosmo:/Volumes/Data2/hestia_muru/muru_frames.json`).

**What was measured, per halo:** the AHF center refined by shrinking spheres
on DM+stars; the disk plane ẑ from the net angular momentum of true stars
(wind excluded) in the 3–15 kpc shell, bulk velocity removed; the stellar
shape-tensor minor axis as an independent cross-check; and the m = 2 bar
amplitude/phase in the face-on disk.

| halo | center shift (kpc) | L̂ vs minor axis | c/a, b/a | A2 max | stars in shell |
|---|---|---|---|---|---|
| G1.1 | 0.021 | 0.22° | 0.56, 0.98 | 0.080 | 2.78M |
| G1.2 | 0.077 | 0.70° | 0.33, 0.95 | 0.149 | 1.99M |
| G2.1 | 0.014 | 0.22° | 0.37, 0.98 | 0.281 | 1.96M |
| G2.2 | 0.038 | 0.07° | 0.48, 0.98 | 0.061 | 1.92M |
| G3.1 | 0.271 | 0.25° | 0.44, 1.00 | 0.079 | 1.08M |
| G3.2 | 0.006 | 0.15° | 0.36, 0.97 | 0.194 | 1.64M |

**Gates.** Center refinement ≤ 0.27 kpc everywhere (gate 1.0 kpc; worst is
G3.1 at ~1.2 softenings — fine). L̂ agrees with the shape-tensor minor axis to
≤ 0.7° in all six — the disk plane is unambiguous, so the Stage 3 observer
placement does not depend on which definition wins.

**Read from the shapes, not asserted:** all six stellar distributions are
oblate disks (b/a ≥ 0.95, c/a 0.33–0.56). The m = 2 amplitude varies a lot:
G2.1 (A2 = 0.28) and G3.2 (0.19) carry real bars; G2.2 (0.06) is nearly
unbarred. So "angle to the bar" is a well-defined view label for some halos
and a weak one for others — Stage 3 records φ_bar per view but analyses must
carry A2 alongside it rather than treating every halo as barred.

⚠ φ_bar is quoted in an arbitrary deterministic in-plane basis (box-x
projected off ẑ); it is a *reference for view angles within a halo*, with no
cross-halo physical meaning.

**Pending cross-check:** Muru's AHF profile files (inertia-tensor axes) were
requested 2026-08-30; when they arrive, compare against `z_hat`/`minor` here.
Our frames stand on their own — the AHF axes are a cross-check, not an input.
