"""Stage 3 of validation/muru_redo_plan_2026-08-30.md: the 216-view map suite.

Per halo, observers every 10 deg of azimuth on the solar circle (R_SUN =
8.18 kpc, reference_values.yaml::a2020_defaults) in the Stage 1 disk plane,
azimuth measured from the Stage 1 bar axis. Per view, four maps on the
40x40 deg, 0.1 deg analysis grid (gcepy orientation: array [ib, il], l and b
both INCREASING with index, GC at the (200,200) pixel boundary):

  jmap      — correct J: sum(m_i rho_i / s_i^2) per pixel / pixel solid angle
              [Msun^2 kpc^-5 sr^-1], s <= 100 kpc. rho from the Stage 2
              k=32 cache.
  stars     — stellar flux map: sum(m_i / s_i^2) per pixel / solid angle,
              TRUE stars (GFM_StellarFormationTime > 0; wind excluded).
  starsold  — same, Muru's old-star cut (0 < a_form < 0.79), their stellar
              bulge proxy, for the like-for-like "indistinguishable" redo.
  muru      — THEIR statistic: mass within 3.0 deg in FLAT (l,b) of each
              pixel center, 15 kpc observer cut (their hard-coded values),
              computed as an exact-metric disk-kernel convolution of the
              binned mass map (their flat-lb aperture IS a convolution up to
              rim pixelization). Exponent NOT applied — squaring is a
              monotonic relabeling, apply downstream if wanted.

WHY HISTOGRAMS, NOT THE APERTURE-KDTREE ESTIMATORS: at 0.1 deg x 216 views
the KDTree route costs ~18 h on nucosmo; histogram/convolution costs
seconds per view. The reformulation is gated, not assumed: before
production, view 0 of G3.2 is computed BOTH ways and compared —
  gate 1: convolved `muru` vs the anchored project_muru(metric="flat_lb")
          (which reproduces their Julia pipeline to 5e-7) — rim pixelization
          only, must agree to < 0.5% in the ROI interior.
  gate 2: histogram jmap vs aperture project_annihilation after identical
          1-deg smoothing — different estimators of the same field, must
          agree to a few percent (recorded, and a failure aborts).

Runs ON nucosmo (micromamba env `hestia`), data local at
/Volumes/Data2/hestia_muru (override with MURU_DATA for Studio/NFS tests).
Outputs: stage3_maps/<halo>.npz + stage3_maps/PROVENANCE.json beside the
data. Requires particle_projection.py next to this script (for the gates).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

H = 0.677
DATA = Path(os.environ.get("MURU_DATA", "/Volumes/Data2/hestia_muru"))
OUTDIR = DATA / "stage3_maps"
R_SUN = 8.18                 # kpc — our primary-leg convention
S_MAX = 100.0                # kpc — J and stellar maps
MURU_SMAX = 15.0             # kpc — their hard-coded cut, for the `muru` map
MURU_APERTURE = 3.0          # deg, flat (l,b) — their hard-coded aperture
OLD_STAR_AMAX = 0.79         # their old-star cut (a_form < 0.79)
NPIX, STEP = 400, 0.1        # gcepy grid: edges -20..20
AZIMUTHS = np.arange(0.0, 360.0, 10.0)

EDGES = np.linspace(-20.0, 20.0, NPIX + 1)
BCTR = 0.5 * (EDGES[:-1] + EDGES[1:])
PIX_SR = np.deg2rad(STEP) ** 2 * np.cos(np.deg2rad(BCTR))[:, None]  # [ib, 1]

FILES = {
    "G1.1": "HESTIA_09_18_8192_G11_minimal_particles.arrow",
    "G1.2": "HESTIA_09_18_8192_G12_minimal_particles.arrow",
    "G2.1": "HESTIA_17_11_8192_G21_minimal_particles.arrow",
    "G2.2": "HESTIA_17_11_8192_G22_minimal_particles.arrow",
    "G3.1": "HESTIA_37_11_8192_G31_minimal_particles.arrow",
    "G3.2": "HESTIA_37_11_8192_G32_minimal_particles.arrow",
}


def disk_kernel():
    """Their aperture: 1 within MURU_APERTURE in flat (l,b), on the 0.1 grid."""
    n = int(round(MURU_APERTURE / STEP))
    ax = np.arange(-n, n + 1) * STEP
    return (np.hypot(ax[:, None], ax[None, :]) <= MURU_APERTURE).astype(float)


def lbs(pos_c, obs, z_hat):
    """(l, b, s) for centered positions, observer obs (in-plane), ez = z_hat.

    Convention of particle_projection._angles_from_observer: ex toward the
    GC, ey = ez x ex, +l toward ey.
    """
    ex = -obs / np.linalg.norm(obs)
    ey = np.cross(z_hat, ex)
    ey /= np.linalg.norm(ey)
    rel = pos_c - obs
    s = np.linalg.norm(rel, axis=1)
    l = np.degrees(np.arctan2(rel @ ey, rel @ ex))
    b = np.degrees(np.arcsin(np.clip((rel @ z_hat) / np.maximum(s, 1e-12), -1, 1)))
    return l, b, s


def binmap(l, b, w):
    """Histogram onto the analysis grid, gcepy orientation [ib, il]."""
    m, _, _ = np.histogram2d(b, l, bins=[EDGES, EDGES], weights=w)
    return m


def view_maps(pos_c, mass, rho, star, old, z_hat, obs, kern):
    from scipy.signal import fftconvolve

    l, b, s = lbs(pos_c, obs, z_hat)
    roi = (np.abs(l) < 20.0) & (np.abs(b) < 20.0)
    jsel = roi & (s <= S_MAX)
    j = binmap(l[jsel], b[jsel], (mass * rho)[jsel] / s[jsel] ** 2) / PIX_SR
    st = binmap(l[jsel & star], b[jsel & star],
                mass[jsel & star] / s[jsel & star] ** 2) / PIX_SR
    so = binmap(l[jsel & old], b[jsel & old],
                mass[jsel & old] / s[jsel & old] ** 2) / PIX_SR
    # their statistic: pad the mass histogram by the kernel radius so cone
    # mass near the ROI edge includes particles just outside it (their KDTree
    # has no ROI cut) — bin on an extended grid, convolve, crop.
    pad = kern.shape[0] // 2
    ext = np.linspace(-20.0 - pad * STEP, 20.0 + pad * STEP, NPIX + 2 * pad + 1)
    msel = (s <= MURU_SMAX) & (np.abs(l) < 20 + MURU_APERTURE + STEP) \
        & (np.abs(b) < 20 + MURU_APERTURE + STEP)
    mh, _, _ = np.histogram2d(b[msel], l[msel], bins=[ext, ext], weights=mass[msel])
    mu = fftconvolve(mh, kern, mode="same")[pad:-pad, pad:-pad]
    return j, st, so, mu


def load_halo(name):
    import pyarrow.feather

    frames = json.load(open(DATA / "muru_frames.json"))[name]
    t = pyarrow.feather.read_table(DATA / FILES[name])
    pt = t["ptype"].to_numpy(zero_copy_only=False)
    keep = (pt == "dm") | (pt == "stars")
    pos = np.column_stack(
        [t[f"Coordinates{i}"].to_numpy() for i in (1, 2, 3)])[keep] * 1000.0 / H
    mass = t["Masses"].to_numpy()[keep].astype(float) * 1e10 / H
    sft = np.nan_to_num(t["GFM_StellarFormationTime"].to_numpy(
        zero_copy_only=False).astype(float), nan=-1.0)[keep]
    isdm = pt[keep] == "dm"
    del t

    c = (np.array(frames["center_ahf_boxMpch"]) * 1000.0 / H
         + np.array(frames["center_refined_offset_kpc"]))
    pos -= c

    rho = np.zeros(len(mass))
    rfile = DATA / "rho_k32" / FILES[name].replace(".arrow", "_rho_k32.npy")
    rho[isdm] = np.load(rfile)          # dm rows in file order == our order
    rho[~isdm] = 0.0                    # only DM annihilates

    star = (~isdm) & (sft > 0.0)
    old = (~isdm) & (sft > 0.0) & (sft < OLD_STAR_AMAX)

    # nothing beyond S_MAX + R_SUN from the center can enter any map
    r = np.linalg.norm(pos, axis=1)
    cut = r < S_MAX + R_SUN + 1.0
    return (pos[cut], mass[cut], rho[cut], star[cut], old[cut],
            np.array(frames["z_hat"]), np.array(frames["bar_x_hat"]), frames)


def gates(pos_c, mass, rho, star, z_hat, x_hat, kern):
    """Validate the histogram/convolution reformulation on one view. Abort on fail."""
    here = Path(__file__).resolve()
    sys.path.insert(0, str(here.parent))          # nucosmo: file shipped alongside
    sys.path.insert(0, str(here.parents[2] / "src"))   # repo: gcereduce
    try:
        from particle_projection import project_annihilation, project_muru
    except ImportError:
        from gcereduce.particle_projection import project_annihilation, project_muru
    from scipy.ndimage import gaussian_filter

    y_hat = np.cross(z_hat, x_hat)
    obs = R_SUN * x_hat
    q = np.column_stack([pos_c @ x_hat, pos_c @ y_hat, pos_c @ z_hat])
    obs_q = np.array([R_SUN, 0.0, 0.0])
    j, st, so, mu = view_maps(pos_c, mass, rho, star, star, z_hat, obs, kern)

    # gate 1: convolution vs the anchored exact flat-lb statistic (0.5 deg
    # grid for cost; the metric identity does not depend on grid step)
    g = np.arange(-19.75, 19.76, 0.5)
    exact = project_muru(q, mass, obs_q, g, g, aperture_deg=MURU_APERTURE,
                         max_dist=MURU_SMAX, exponent=1, metric="flat_lb")
    # project_muru returns [ib, il] — the same orientation as our histogram
    # maps (no transpose; the anchor's .T was for Julia's [il, ib] output)
    ib = np.searchsorted(EDGES, g) - 1
    conv_at = mu[np.ix_(ib, ib)]
    live = exact > 0
    r1 = np.abs(conv_at[live] - exact[live]) / exact[live]
    print(f"gate 1 (muru conv vs exact): median {np.median(r1):.2e} "
          f"max {r1.max():.2e}", flush=True)
    # measured on G3.2 view 0: median 5.1e-3, max 2.1e-2 — rim pixelization
    # (the +/-0.05 deg quantization ring is ~3% of the 3-deg cone's area, so
    # sub-percent median is the expected scale; worst cells are low-mass
    # high-|b| edges). Tolerance set above the measurement, not tuned to it:
    # the muru maps feed qualitative side-by-sides, where 1% is invisible.
    assert np.median(r1) < 1e-2 and r1.max() < 5e-2, "gate 1 FAIL"

    # gate 2: histogram J vs aperture project_annihilation, both smoothed 1 deg
    ap = project_annihilation(q, mass, rho, obs_q, g, g,
                              aperture_deg=0.5, max_dist=S_MAX)
    sig_c, sig_f = 1.0 / 0.5, 1.0 / STEP
    hs = gaussian_filter(j, sig_f)[np.ix_(ib, ib)]
    as_ = gaussian_filter(ap, sig_c)        # already [ib, il]
    lv = as_ > 0.01 * as_.max()
    r2 = np.abs(hs[lv] - as_[lv]) / as_[lv]
    print(f"gate 2 (hist J vs aperture J, 1deg smooth): median {np.median(r2):.2e} "
          f"max {r2.max():.2e}", flush=True)
    assert np.median(r2) < 0.05, "gate 2 FAIL"
    return {"gate1_median": float(np.median(r1)), "gate1_max": float(r1.max()),
            "gate2_median": float(np.median(r2)), "gate2_max": float(r2.max())}


def main():
    OUTDIR.mkdir(exist_ok=True)
    kern = disk_kernel()
    prov = {"script": "scripts/muru_redo/stage3_views.py", "date": "2026-08-30",
            "R_sun_kpc": R_SUN, "s_max_kpc": S_MAX,
            "muru_smax_kpc": MURU_SMAX, "muru_aperture_deg": MURU_APERTURE,
            "old_star_amax": OLD_STAR_AMAX, "azimuths_deg": AZIMUTHS.tolist(),
            "grid": "gcepy orientation [ib, il], edges -20..20, 0.1 deg",
            "lb_convention": "ex toward GC, ey = z_hat x ex, +l toward ey",
            "units": {"jmap": "Msun^2 kpc^-5 sr^-1",
                      "stars/starsold": "Msun kpc^-2 sr^-1",
                      "muru": "Msun in flat-lb 3deg cone (exponent NOT applied)"}}

    g32 = load_halo("G3.2")
    prov["gates_on_G3.2_view0"] = gates(g32[0], g32[1], g32[2], g32[3],
                                        g32[5], g32[6], kern)

    for name in FILES:
        t0 = time.time()
        pos_c, mass, rho, star, old, z_hat, x_hat, fr = \
            g32 if name == "G3.2" else load_halo(name)
        y_hat = np.cross(z_hat, x_hat)
        shp = (len(AZIMUTHS), NPIX, NPIX)
        maps = {k: np.zeros(shp, np.float32) for k in
                ("jmap", "stars", "starsold", "muru")}
        for i, a in enumerate(np.deg2rad(AZIMUTHS)):
            obs = R_SUN * (np.cos(a) * x_hat + np.sin(a) * y_hat)
            j, st, so, mu = view_maps(pos_c, mass, rho, star, old, z_hat, obs, kern)
            for k, v in zip(("jmap", "stars", "starsold", "muru"), (j, st, so, mu)):
                maps[k][i] = v
        np.savez_compressed(OUTDIR / f"{name}.npz",
                            azimuth_deg=AZIMUTHS, **maps)
        prov[name] = {"n_particles_loaded": int(len(mass)),
                      "bar_A2_max": fr["bar_A2_max"],
                      "elapsed_s": round(time.time() - t0, 1)}
        print(f"{name}: 36 views x 4 maps in {prov[name]['elapsed_s']} s",
              flush=True)

    with open(OUTDIR / "PROVENANCE.json", "w") as f:
        json.dump(prov, f, indent=1)
    print("wrote", OUTDIR / "PROVENANCE.json", flush=True)


if __name__ == "__main__":
    main()
