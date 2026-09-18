"""Stage 2 of validation/muru_redo_plan_2026-08-30.md: kNN density estimation.

Part A (--converge): k-convergence on one halo (G3.2, the anchored one) BEFORE
production. Two measurements, k in {16, 32, 64}:
  1. per-particle rho_k / rho_32 in radial bins — estimator systematics vs
     radius, where softening and particle noise live;
  2. the quantity we actually consume: correct J-maps (project_annihilation)
     from one solar-circle view, smoothed to the Stage 4 working resolution,
     max/median relative differences.

Part B (--produce): rho at the adopted k for all six halos, cached float32
beside the data (rho_k32/*.npy, ~40 MB total). rho is translation-invariant,
so no centering enters — positions are only scaled to physical kpc; units
Msun/kpc^3.

Runs on the Studio against the NFS mount (Kev, 2026-08-30):
/private/nfs/Data2/hestia_muru. Reads frames from muru_frames.json (Stage 1).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from particle_projection import knn_density, project_annihilation  # noqa: E402

H = 0.677
DATA = Path("/private/nfs/Data2/hestia_muru")
RHODIR = DATA / "rho_k32"
R_SUN = 8.18            # kpc, reference_values.yaml::a2020_defaults (primary leg)
S_MAX = 100.0           # kpc
GRID = np.arange(-20.0, 20.0 + 1e-9, 0.5)
APERTURE = 0.5          # deg — convergence-check beams
SMOOTH_SIGMA_DEG = 1.0  # Gaussian sigma; Stage 4 reports at 2-3 deg FWHM-ish scales

FILES = {
    "G1.1": "HESTIA_09_18_8192_G11_minimal_particles.arrow",
    "G1.2": "HESTIA_09_18_8192_G12_minimal_particles.arrow",
    "G2.1": "HESTIA_17_11_8192_G21_minimal_particles.arrow",
    "G2.2": "HESTIA_17_11_8192_G22_minimal_particles.arrow",
    "G3.1": "HESTIA_37_11_8192_G31_minimal_particles.arrow",
    "G3.2": "HESTIA_37_11_8192_G32_minimal_particles.arrow",
}


def load_dm(name):
    """DM positions (physical kpc, box-origin) and masses (Msun)."""
    import pyarrow.feather

    t = pyarrow.feather.read_table(DATA / FILES[name])
    dm = t["ptype"].to_numpy(zero_copy_only=False) == "dm"
    pos = np.column_stack(
        [t[f"Coordinates{i}"].to_numpy() for i in (1, 2, 3)])[dm] * 1000.0 / H
    mass = t["Masses"].to_numpy()[dm].astype(float) * 1e10 / H
    return pos, mass


def centered(name, pos):
    """Positions relative to the Stage 1 refined center, physical kpc."""
    frames = json.load(open(DATA / "muru_frames.json"))
    f = frames[name]
    c = (np.array(f["center_ahf_boxMpch"]) * 1000.0 / H
         + np.array(f["center_refined_offset_kpc"]))
    return pos - c, f


def view_jmap(pos_c, mass, rho, frame):
    """Correct J-map from one view: observer at R_SUN along the bar axis."""
    z = np.array(frame["z_hat"])
    x = np.array(frame["bar_x_hat"])
    y = np.cross(z, x)
    q = pos_c @ np.column_stack([x, y, z])          # disk frame
    obs = np.array([R_SUN, 0.0, 0.0])
    return project_annihilation(q, mass, rho, obs, GRID, GRID,
                                aperture_deg=APERTURE, max_dist=S_MAX)


def smooth(m):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(m, SMOOTH_SIGMA_DEG / (GRID[1] - GRID[0]))


def converge():
    name = "G3.2"
    pos, mass = load_dm(name)
    pos_c, frame = centered(name, pos)
    r = np.linalg.norm(pos_c, axis=1)

    rhos, jmaps = {}, {}
    for k in (16, 32, 64):
        t0 = time.time()
        rhos[k] = knn_density(pos_c, mass, k=k)
        jmaps[k] = view_jmap(pos_c, mass, rhos[k], frame)
        print(f"k={k}: rho+jmap in {time.time()-t0:.0f} s", flush=True)

    print("\nper-particle rho_k/rho_32, median [16-84%], by radius:")
    bins = [(0.22, 1.0), (1.0, 3.0), (3.0, 10.0), (10.0, 30.0)]
    for lo, hi in bins:
        s = (r > lo) & (r < hi)
        for k in (16, 64):
            q = rhos[k][s] / rhos[32][s]
            p16, p50, p84 = np.percentile(q, [16, 50, 84])
            print(f"  {lo:5.2f}-{hi:5.2f} kpc  k={k:2d}: "
                  f"{p50:.4f} [{p16:.4f}, {p84:.4f}]  (n={s.sum():,})")

    print(f"\nsmoothed J-map (sigma={SMOOTH_SIGMA_DEG} deg) rel diff vs k=32, "
          f"pixels > 1% of max:")
    j32 = smooth(jmaps[32])
    live = j32 > 0.01 * j32.max()
    out = {}
    for k in (16, 64):
        d = np.abs(smooth(jmaps[k])[live] - j32[live]) / j32[live]
        out[k] = (float(np.median(d)), float(d.max()))
        print(f"  k={k:2d}: median {out[k][0]:.4f}  max {out[k][1]:.4f}  "
              f"({live.sum()} px)")
    return out


def produce():
    RHODIR.mkdir(exist_ok=True)
    prov = {"script": "scripts/muru_redo/stage2_density.py", "k": 32,
            "date": "2026-08-30", "units": "Msun/kpc^3 (physical, h removed)",
            "note": "rho for DM particles in file order (dm rows only)"}
    for name, fname in FILES.items():
        t0 = time.time()
        pos, mass = load_dm(name)
        rho = knn_density(pos, mass, k=32)
        out = RHODIR / fname.replace(".arrow", "_rho_k32.npy")
        np.save(out, rho.astype(np.float32))
        prov[name] = {"file": out.name, "n": len(rho),
                      "rho_max": float(rho.max()), "rho_median": float(np.median(rho))}
        print(f"{name}: {len(rho):,} rho in {time.time()-t0:.0f} s -> {out.name}",
              flush=True)
    with open(RHODIR / "PROVENANCE.json", "w") as f:
        json.dump(prov, f, indent=1)
    print("wrote", RHODIR / "PROVENANCE.json")


if __name__ == "__main__":
    if "--produce" in sys.argv:
        produce()
    else:
        converge()
