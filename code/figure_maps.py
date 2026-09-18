"""Generate the maps for the Muru figure reproduction (runs ON nucosmo).

Their Figs. 1-2 geometry (arXiv:2508.06314 Appendix A + their code, anchored
2026-08-30): observer in the disk plane, 8.0 kpc from the center, along the
galaxy's MAJOR AXIS; 3-deg flat-(l,b) aperture; 15 kpc distance cut; maps
normalized to max. Per halo, four maps on the 0.1-deg grid:

  mu_dm     — their linear DM panel: cone mass, dm particles (exponent 1;
              their "quadratic" panel is this map squared — a monotonic
              relabeling, drawn downstream, not recomputed).
  mu_stars  — their stellar panel: same statistic on their old-star cut
              (0 < a_form < 0.79).
  jmap      — the corrected annihilation map from the SAME viewpoint
              (sum m rho / s^2 per pixel, s <= 100 kpc, Stage 2 rho).
  jmap_muru_geom — corrected J with THEIR 15 kpc truncation, for a
              like-for-like corrected panel.

FRAMES: **AHF, as of 2026-09-15.** Muru delivered the AHF profile files, so
the viewing frame is now READ from them (`ahf_frames.py`, Ea = observer axis,
Ec = disk normal, at the converged bin nearest 10 kpc/h) rather than
reconstructed. `their_frame()` below is kept as the fallback and is still
evaluated every run, so the log carries the angle between the two.

Set MURU_FRAME_SOURCE=particles to force the old reconstruction.

Everything else is their pipeline, anchored to 5e-7. Output: figure_maps.npz
beside the data.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
DATA = Path(os.environ.get("MURU_DATA", "/Volumes/Data2/hestia_muru"))
sys.path.insert(0, str(DATA))

import stage3_views as s3  # noqa: E402  (same dir on nucosmo)

D_GC = 8.0        # kpc — THEIR observer distance (not our 8.18)
MURU_SMAX = 15.0  # kpc — their cut


def their_frame(pos_c, mass):
    """THEIR frame, per the published Appendix A steps 1-2 (PRL 135, 161005):
    inertia tensor from particles <= 10/h kpc = 14.8 kpc of the center (their
    files ARE the AHF particle lists, all species), minor axis = disk normal,
    major axis = observer direction. Moment-of-inertia and shape tensors
    share eigenvectors (I = tr(S)delta - S), so the shape tensor's smallest/
    largest eigenvectors are the minor/major axes. Sign of each axis is
    ambiguous (+/-); the axis-ratio statistics are mirror-invariant to
    <= 0.004 (Stage 4), so the comparison is insensitive to it."""
    r = np.linalg.norm(pos_c, axis=1)
    s = r < 10.0 / 0.677
    w = mass[s]
    q = pos_c[s]
    S = np.einsum("i,ij,ik->jk", w, q, q) / w.sum()
    ev, evec = np.linalg.eigh(S)
    v_n, v_major = evec[:, 0], evec[:, 2]
    v_major = v_major - (v_major @ v_n) * v_n       # exactly orthogonal
    return v_n, v_major / np.linalg.norm(v_major)


def ahf_frames():
    """Delivered AHF axes, keyed by halo. Written by `ahf_frames.py`; copied
    next to the data so this script stays runnable on nucosmo."""
    # MURU_AHF_FRAMES selects the frames file (default: observer on the AHF
    # major axis). "ahf_frames_bestaz.json" is option (b), 2026-09-16: the
    # in-plane azimuth chosen per halo, within the b/a ~ 0.99 degeneracy, to
    # match the published corner tables (ahf_azimuth_scan.py, 10 deg steps).
    fname = os.environ.get("MURU_AHF_FRAMES", "ahf_frames.json")
    for cand in (DATA / fname,
                 Path(__file__).resolve().parents[1] / "figures" / fname):
        if cand.exists():
            return json.load(open(cand))
    raise FileNotFoundError("ahf_frames.json not found; run analysis/code/ahf_frames.py")


def axis_angle(u, v):
    c = abs(float(np.dot(u, v))) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(min(1.0, c))))


def load_gas(name, c_ahf_kpc):
    """Gas rows only (dropped by stage3's loader), AHF-centered physical kpc.
    Needed because their AHF inertia tensor uses ALL bound particles."""
    import pyarrow.feather

    t = pyarrow.feather.read_table(DATA / s3.FILES[name])
    g = t["ptype"].to_numpy(zero_copy_only=False) == "gas"
    pos = np.column_stack(
        [t[f"Coordinates{i}"].to_numpy() for i in (1, 2, 3)])[g] * 1000.0 / 0.677
    return pos - c_ahf_kpc, t["Masses"].to_numpy()[g].astype(float) * 1e10 / 0.677


def main():
    use_ahf = os.environ.get("MURU_FRAME_SOURCE", "ahf") == "ahf"
    AHF = ahf_frames() if use_ahf else None
    out, meta = {}, {"D_GC_kpc": D_GC,
                     "axis": ("AHF profile axes as delivered by Muru "
                              "2026-09-15: Ea = observer, Ec = normal, at the "
                              "converged bin nearest 10 kpc/h" if use_ahf else
                              "reconstructed: all-particle (incl. gas) shape "
                              "tensor <= 10/h kpc about the AHF center"),
                     "frame_source": "ahf" if use_ahf else "particles",
                     "frames_file": os.environ.get("MURU_AHF_FRAMES", "ahf_frames.json") if use_ahf else None,
                     "date": "2026-09-15"}
    kern = s3.disk_kernel()
    frames = json.load(open(DATA / "muru_frames.json"))
    for name in s3.FILES:
        t0 = time.time()
        pos_c, mass, rho, star, old, z_hat, x_hat, fr = s3.load_halo(name)
        off = np.array(frames[name]["center_refined_offset_kpc"])
        pos_c = pos_c + off                       # back to the AHF center
        c_ahf = np.array(frames[name]["center_ahf_boxMpch"]) * 1000.0 / 0.677
        gpos, gmass = load_gas(name, c_ahf)
        tpos = np.vstack([pos_c, gpos])
        tmass = np.concatenate([mass, gmass])
        v_n_rec, maj_rec = their_frame(tpos, tmass)
        del tpos, tmass, gpos, gmass
        if use_ahf:
            a = AHF[name]
            v_n, maj = np.array(a["v_normal"]), np.array(a["v_major"])
            d_maj = axis_angle(maj, maj_rec)
            d_nrm = axis_angle(v_n, v_n_rec)
            print(f"  {name}: AHF frame (r={a['r_kpch']:.2f} kpc/h, az={a.get('azimuth_deg', 0)}, "
                  f"b/a={a['b_over_a']:.3f}); reconstruction differs by "
                  f"{d_maj:.1f} deg (major), {d_nrm:.1f} deg (normal)")
        else:
            v_n, maj, d_maj, d_nrm = v_n_rec, maj_rec, 0.0, 0.0
        obs = D_GC * maj
        l, b, s = s3.lbs(pos_c, obs, v_n)

        pad = kern.shape[0] // 2
        ext = np.linspace(-20 - pad * s3.STEP, 20 + pad * s3.STEP,
                          s3.NPIX + 2 * pad + 1)

        def conemass(sel):
            from scipy.signal import fftconvolve
            m = (s <= MURU_SMAX) & sel & \
                (np.abs(l) < 20 + 3.0 + s3.STEP) & (np.abs(b) < 20 + 3.0 + s3.STEP)
            h, _, _ = np.histogram2d(b[m], l[m], bins=[ext, ext], weights=mass[m])
            return fftconvolve(h, kern, mode="same")[pad:-pad, pad:-pad]

        dm_sel = rho > 0                            # rho only set for dm rows
        roi = (np.abs(l) < 20) & (np.abs(b) < 20)

        def jmap(smax):
            sel = roi & (s <= smax)
            return s3.binmap(l[sel], b[sel],
                             (mass * rho)[sel] / s[sel] ** 2) / s3.PIX_SR

        def colmap(which, smax):
            """Column density int rho ds = sum(m_i/s_i^2)/dOmega per pixel.

            🚨 THE CORRECT STELLAR COMPARISON MAP (Jason Kumar, 2026-08-31).
            A stellar source at distance s dilutes as 1/s^2 exactly as
            annihilation emission does, and the beam volume grows as s^2, so
            the two cancel and the observable morphology of an old stellar
            population is the COLUMN DENSITY int rho_* ds -- not the
            s^2-weighted cone mass int rho_* s^2 ds that Ref. [1] plots. The
            cone-mass version is kept as `mu_stars` for reproducing their
            panels; this is what the corrected annihilation map must be
            compared against.
            """
            sel = roi & which & (s <= smax)
            return s3.binmap(l[sel], b[sel], mass[sel] / s[sel] ** 2) / s3.PIX_SR

        out[f"{name}_mu_dm"] = conemass(dm_sel).astype(np.float32)
        out[f"{name}_mu_stars"] = conemass(old).astype(np.float32)
        out[f"{name}_stars_col"] = colmap(old, MURU_SMAX).astype(np.float32)
        out[f"{name}_jmap"] = jmap(100.0).astype(np.float32)
        out[f"{name}_jmap_muru_geom"] = jmap(MURU_SMAX).astype(np.float32)
        meta[name] = {"major_axis": [float(v) for v in maj],
                      "normal_axis": [float(v) for v in v_n],
                      "major_axis_reconstructed": [float(v) for v in maj_rec],
                      "normal_axis_reconstructed": [float(v) for v in v_n_rec],
                      "d_major_ahf_vs_reconstructed_deg": d_maj,
                      "d_normal_ahf_vs_reconstructed_deg": d_nrm,
                      "bar_A2_max": fr["bar_A2_max"],
                      "angle_to_bar_deg": float(np.degrees(np.arccos(
                          np.clip(abs(maj @ np.array(fr["bar_x_hat"])), 0, 1)))),
                      "elapsed_s": round(time.time() - t0, 1)}
        print(f"{name}: done [{meta[name]['elapsed_s']} s], "
              f"axis-to-bar {meta[name]['angle_to_bar_deg']:.1f} deg", flush=True)

    np.savez_compressed(DATA / "figure_maps.npz", **out)
    with open(DATA / "figure_maps_meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print("wrote", DATA / "figure_maps.npz", flush=True)


if __name__ == "__main__":
    main()
