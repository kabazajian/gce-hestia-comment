"""Stage 1 of validation/muru_redo_plan_2026-08-30.md: per-halo Galactic frames.

For each of the six HESTIA halos:

1. CENTER — AHF center (Muru email 2026-08-30), refined by shrinking spheres
   on DM+stars. Gate: the refinement must move the center by less than
   GATE_CENTER_KPC (a few softenings; paper softening 0.22 kpc) or the halo is
   flagged rather than silently accepted.
2. DISK PLANE — z_hat = direction of the net angular momentum of TRUE stars
   (GFM_StellarFormationTime > 0; <= 0 is stellar wind) in the 3-15 kpc
   physical shell, bulk velocity removed. Cross-checked against the stellar
   inertia (shape) tensor minor axis; misalignment > FLAG_MISALIGN_DEG flags
   the halo.
3. BAR — m=2 Fourier mode of the face-on stellar surface density in
   cylindrical bins; phi_bar from the mass-weighted m=2 phase over the bar
   region, A2 profile recorded. The bar axis is the in-plane x_hat reference
   every Stage 3 view angle is quoted against.

Runs on nucosmo (micromamba env `hestia`, data in /Volumes/Data2/hestia_muru).
Output: muru_frames.json next to the data (copied back into the repo by the
launcher). Pure numpy + pyarrow; no repo imports, so it can run standalone.

Conventions: physical kpc via (x_boxMpc/h - c) * 1000 / h at z=0; masses in
Msun via m * 1e10 / h; velocities km/s (z=0, sqrt(a)=1).
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

H = 0.677  # HESTIA cosmology, Libeskind et al. 2018 (Muru email 2026-08-30)
DATA_DIR = "/Volumes/Data2/hestia_muru"
OUT_JSON = f"{DATA_DIR}/muru_frames.json"

# AHF centers, absolute box Mpc/h (Muru email 2026-08-30)
CENTERS = {
    "G1.1": ("HESTIA_09_18_8192_G11_minimal_particles.arrow", (47.3092422, 48.8026602, 50.0027344)),
    "G1.2": ("HESTIA_09_18_8192_G12_minimal_particles.arrow", (46.7936523, 49.0552930, 49.8811797)),
    "G2.1": ("HESTIA_17_11_8192_G21_minimal_particles.arrow", (48.8151992, 46.7045156, 53.6050664)),
    "G2.2": ("HESTIA_17_11_8192_G22_minimal_particles.arrow", (48.7157578, 47.0623281, 53.3371875)),
    "G3.1": ("HESTIA_37_11_8192_G31_minimal_particles.arrow", (46.3911250, 50.7461367, 47.9353516)),
    "G3.2": ("HESTIA_37_11_8192_G32_minimal_particles.arrow", (46.7533242, 50.3214297, 47.7923984)),
}

GATE_CENTER_KPC = 1.0      # ~4-5 softenings; larger => flag, do not proceed silently
FLAG_MISALIGN_DEG = 10.0   # L_hat vs inertia minor axis
SHELL_LO, SHELL_HI = 3.0, 15.0   # kpc, disk-plane shell
BAR_RMAX, BAR_ZMAX = 6.0, 3.0    # kpc, bar search region (face-on cut applied)


def shrinking_spheres(pos_kpc, mass, r0=20.0, shrink=0.9, rmin=0.7, nmin=1000):
    """Refine the center: iterate COM of particles within r, shrinking r."""
    c = np.zeros(3)
    r = r0
    d = np.linalg.norm(pos_kpc - c, axis=1)
    sel = d < r
    while True:
        w = mass[sel]
        c_new = (pos_kpc[sel] * w[:, None]).sum(axis=0) / w.sum()
        r *= shrink
        d = np.linalg.norm(pos_kpc - c_new, axis=1)
        sel = d < r
        c = c_new
        if r < rmin or sel.sum() < nmin:
            return c


def frame_for(name, fname, center_boxMpch):
    import read_particles

    t0 = time.time()
    t = read_particles.read_table(f"{DATA_DIR}/{fname}")
    pt = t["ptype"].to_numpy(zero_copy_only=False)
    pos = np.column_stack([t[f"Coordinates{i}"].to_numpy() for i in (1, 2, 3)])
    vel = np.column_stack([t[f"Velocities{i}"].to_numpy() for i in (1, 2, 3)]).astype(float)
    mass = t["Masses"].to_numpy().astype(float) * 1e10 / H          # Msun
    sft = np.nan_to_num(
        t["GFM_StellarFormationTime"].to_numpy(zero_copy_only=False).astype(float), nan=-1.0)
    del t

    pos_kpc = (pos - np.asarray(center_boxMpch)) * 1000.0 / H       # physical kpc
    del pos
    core = pt != "gas"                                              # DM + stars (+bh)
    star = (pt == "stars") & (sft > 0.0)

    # --- 1. center ---------------------------------------------------------
    c_ref = shrinking_spheres(pos_kpc[core], mass[core])
    dc = float(np.linalg.norm(c_ref))
    pos_kpc -= c_ref                                                # refined center = origin

    # --- 2. disk plane -----------------------------------------------------
    r = np.linalg.norm(pos_kpc, axis=1)
    shell = star & (r > SHELL_LO) & (r < SHELL_HI)
    w = mass[shell]
    vbulk = (vel[shell] * w[:, None]).sum(axis=0) / w.sum()
    L = np.cross(pos_kpc[shell], vel[shell] - vbulk)
    Lg = (L * w[:, None]).sum(axis=0)
    z_hat = Lg / np.linalg.norm(Lg)

    # inertia (shape) tensor of the same stars; minor axis = smallest eigval
    q = pos_kpc[shell]
    S = np.einsum("i,ij,ik->jk", w, q, q) / w.sum()
    evals, evecs = np.linalg.eigh(S)
    minor = evecs[:, 0]
    mis = float(np.degrees(np.arccos(min(1.0, abs(minor @ z_hat)))))
    axis_ratios = [float(np.sqrt(evals[0] / evals[2])), float(np.sqrt(evals[1] / evals[2]))]

    # --- 3. bar ------------------------------------------------------------
    # in-plane basis (deterministic): x0 = box-x projected off z_hat
    x0 = np.array([1.0, 0.0, 0.0])
    x0 = x0 - (x0 @ z_hat) * z_hat
    x0 /= np.linalg.norm(x0)
    y0 = np.cross(z_hat, x0)
    px, py, pz = pos_kpc @ x0, pos_kpc @ y0, pos_kpc @ z_hat
    rcyl = np.hypot(px, py)
    disk = star & (rcyl < BAR_RMAX) & (np.abs(pz) < BAR_ZMAX)
    phi = np.arctan2(py[disk], px[disk])
    wm = mass[disk]
    edges = np.arange(0.0, BAR_RMAX + 1e-9, 0.5)
    a2_prof, phi_prof = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (rcyl[disk] >= lo) & (rcyl[disk] < hi)
        if s.sum() < 100:
            a2_prof.append(0.0); phi_prof.append(np.nan); continue
        z2 = (wm[s] * np.exp(2j * phi[s])).sum() / wm[s].sum()
        a2_prof.append(float(np.abs(z2)))
        phi_prof.append(float(np.degrees(0.5 * np.angle(z2))))
    a2_prof = np.array(a2_prof)
    a2_max = float(a2_prof.max())
    # bar phase: mass-weighted m=2 over the region where A2 > half its peak,
    # restricted to r < 5 kpc (a bar is an inner feature; spiral arms are not)
    inner = rcyl[disk] < 5.0
    z2all = (wm[inner] * np.exp(2j * phi[inner])).sum() / wm[inner].sum()
    phi_bar = float(np.degrees(0.5 * np.angle(z2all)))
    bar_x = np.cos(np.radians(phi_bar)) * x0 + np.sin(np.radians(phi_bar)) * y0

    flags = []
    if dc > GATE_CENTER_KPC:
        flags.append(f"center moved {dc:.2f} kpc > {GATE_CENTER_KPC} gate")
    if mis > FLAG_MISALIGN_DEG:
        flags.append(f"L/inertia misalignment {mis:.1f} deg > {FLAG_MISALIGN_DEG}")

    out = {
        "file": fname,
        "center_ahf_boxMpch": list(center_boxMpch),
        "center_refined_offset_kpc": [float(v) for v in c_ref],
        "center_shift_kpc": dc,
        "z_hat": [float(v) for v in z_hat],
        "L_star_msun_kpc_kms": float(np.linalg.norm(Lg)),
        "n_stars_shell": int(shell.sum()),
        "vbulk_kms": [float(v) for v in vbulk],
        "inertia_minor_axis": [float(v) for v in minor],
        "inertia_axis_ratios_ca_ba": axis_ratios,
        "misalign_L_vs_minor_deg": mis,
        "bar_x_hat": [float(v) for v in bar_x],
        "bar_phi_deg_in_x0y0": phi_bar,
        "bar_A2_profile_r0.5kpc_bins": [float(v) for v in a2_prof],
        "bar_A2_max": a2_max,
        "flags": flags,
        "elapsed_s": round(time.time() - t0, 1),
    }
    print(f"{name}: center shift {dc:.3f} kpc | misalign {mis:.2f} deg | "
          f"A2max {a2_max:.3f} | phi_bar {phi_bar:+.1f} deg | flags {flags or 'none'} "
          f"[{out['elapsed_s']} s]", flush=True)
    return out


def main():
    sys.path.insert(0, DATA_DIR)
    results = {"_provenance": {
        "script": "scripts/muru_redo/stage1_frames.py",
        "date": "2026-08-30",
        "h": H,
        "shell_kpc": [SHELL_LO, SHELL_HI],
        "gate_center_kpc": GATE_CENTER_KPC,
        "star_def": "ptype==stars and GFM_StellarFormationTime>0 (wind excluded)",
    }}
    for name, (fname, c) in CENTERS.items():
        results[name] = frame_for(name, fname, c)
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
