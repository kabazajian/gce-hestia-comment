"""Stage 3a anchor comparison: our `project_muru` vs Muru et al.'s own code.

Run AFTER anchor_muru_g32.jl. Three layers, each falsifiable:

1. TRANSFORM: our implementation of their (l, b, d) transform must match their
   `add_galacticcoordinates!` output per particle (julia_lbd_head.tsv).
2. MAP: our `project_muru(metric="flat_lb", exponent=1)` on the pre-rotated
   frame must match their `calc_density_on_grid` map cell by cell.
3. DISCRIMINATION: the match must FAIL when the ℓ axis is mirrored and when
   the aperture metric is changed to "flat_lcosb" — a comparison that cannot
   fail proves nothing (project rule; and this is coordinate-trap country,
   CONVENTIONS.md counts seven).

Geometry constants are duplicated from anchor_muru_g32.jl — keep in sync.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from particle_projection import project_muru  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
DATA = Path("/Users/aba/work/hestia_muru/HESTIA_37_11_8192_G32_minimal_particles.arrow")

H = 0.677
P_C = np.array([46.7533242, 50.3214297, 47.7923984]) * 1e3          # kpc/h
V_N = np.array([0.2, 0.3, 0.9]); V_N = V_N / np.linalg.norm(V_N)
_vd = np.array([1.0, 0.1, -0.05]); _vd = _vd - (_vd @ V_N) * V_N
V_D = _vd / np.linalg.norm(_vd)
ANGLE, D_GC, D_MAX, APERTURE = 0.0, 8.0 * H, 15.0 * H, 3.0
GRID = np.arange(-20.0, 20.0 + 1e-9, 0.1)                            # their -20:0.1:20


def their_transform(pos_kpch, p_obs):
    """Literal numpy port of GalacticCoordinates.convert_xyz_to_lb."""
    e_gc = P_C - p_obs
    e_gc = e_gc / np.linalg.norm(e_gc)
    e_n = V_N
    e_y = np.cross(e_gc, e_n)
    e_y = e_y / np.linalg.norm(e_y)
    R = np.vstack([e_gc, e_y, e_n])
    rel = (pos_kpch - p_obs) @ R.T
    d = np.linalg.norm(rel, axis=1)
    l = np.degrees(np.arctan2(rel[:, 1], rel[:, 0]))
    b = np.degrees(np.arcsin(np.clip(rel[:, 2] / np.maximum(d, 1e-12), -1, 1)))
    return l, b, d, R


def main():
    import pyarrow.feather  # Arrow IPC file

    t = pyarrow.feather.read_table(DATA)
    pt = t["ptype"].to_numpy(zero_copy_only=False)
    dm = pt == "dm"
    pos = np.column_stack([t[f"Coordinates{i}"].to_numpy() for i in (1, 2, 3)])[dm] * 1e3
    mass = t["Masses"].to_numpy()[dm].astype(float)
    pids = t["ParticleIDs"].to_numpy()[dm]
    assert dm.sum() == 4_472_132, f"dm count {dm.sum()}"

    # observer: generate_sun_position with alpha=0 -> p_c + D_GC * v_d
    p_obs = P_C + D_GC * V_D
    j_obs = np.loadtxt(OUT / "julia_observer.tsv")
    assert np.allclose(p_obs, j_obs, rtol=0, atol=1e-6), f"observer differs: {p_obs} vs {j_obs}"

    l, b, d, R = their_transform(pos, p_obs)

    # ---- layer 1: per-particle transform vs their output --------------------
    ref = np.loadtxt(OUT / "julia_lbd_head.tsv")
    idx = {p: i for i, p in enumerate(pids)}
    rows = np.array([idx[int(p)] for p in ref[:, 0]])
    for k, (name, ours) in enumerate((("l", l), ("b", b), ("d", d)), start=1):
        diff = np.abs(ours[rows] - ref[:, k]).max()
        print(f"transform {name}: max |ours - theirs| = {diff:.3e}")
        assert diff < 1e-4, f"{name} transform mismatch"

    # ---- layer 2: the map through OUR library -------------------------------
    # Pre-rotate into their frame: GC at origin, observer at (-D_GC, 0, 0).
    q = (pos - P_C) @ R.T
    obs_q = np.array([-D_GC, 0.0, 0.0])
    ours = project_muru(q, mass, obs_q, GRID, GRID,
                        aperture_deg=APERTURE, max_dist=D_MAX,
                        exponent=1, metric="flat_lb")
    theirs = np.loadtxt(OUT / "julia_losmass.tsv")   # [i, j] = (l_i, b_j)

    # theirs indexes (l, b); ours from meshgrid is [b, a] = (b_j, l_a), so a
    # plain transpose aligns them. NO l-flip: the particles were pre-rotated
    # with THEIR R matrix, so both codes share axes by construction. (The
    # library's own-frame mirror — our e_y = e_n x e_gc vs their e_gc x e_n —
    # applies only when each code builds its frame from raw coordinates; a
    # first version of this script applied it here too and failed at 26%.)
    ours_lb = ours.T
    live = theirs > 0
    rel = np.abs(ours_lb[live] - theirs[live]) / theirs[live]
    print(f"map: {live.sum()} live cells; max rel diff = {rel.max():.3e}; "
          f"median = {np.median(rel):.3e}")
    assert rel.max() < 1e-5, "map mismatch"     # residual = their Float32 accumulation

    # ---- layer 3: the comparison can fail -----------------------------------
    flipped = ours[:, ::-1].T
    bad = np.abs(flipped[live] - theirs[live]) / theirs[live]
    print(f"discrimination (l-flipped): max rel diff = {bad.max():.3e}")
    assert bad.max() > 0.1, "l-flip does not discriminate — comparison is vacuous"

    lcosb = project_muru(q, mass, obs_q, GRID, GRID,
                         aperture_deg=APERTURE, max_dist=D_MAX,
                         exponent=1, metric="flat_lcosb").T
    m = np.abs(lcosb[live] - theirs[live]) / theirs[live]
    print(f"discrimination (flat_lcosb metric): max rel diff = {m.max():.3e} "
          f"(their flat-lb aperture is not a constant angular aperture)")
    assert m.max() > 1e-3, "metric choice does not discriminate"

    print("ANCHOR PASS: project_muru(metric='flat_lb') reproduces their pipeline")


if __name__ == "__main__":
    main()
