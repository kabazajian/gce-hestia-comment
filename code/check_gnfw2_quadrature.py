"""Pixel-quadrature check on the analytic gNFW^2 template (2026-09-08).

WHY. `stage4_templates.py` builds the gNFW^2 comparison map by evaluating
`dm_templates.j_map` at the CENTRES of the 0.1 deg analysis pixels. The LOS
integral itself is fine -- it is a genuine int rho^2 ds on sinh-substituted
nodes, converged to 1e-6 in n_s -- but the SKY-PLANE quadrature is a midpoint
rule on a gamma = 1.2 cusp, which is the failure mode `j_factor_roi`'s own
docstring warns about for Cartesian grids. It matters here because ~14% of the
flux inside 3 deg lies inside 0.1 deg: the central pixel comes out 1.7x low.

WHAT THIS DOES. The spherical J field depends only on the great-circle angle
psi, so J(psi) is computed once on a fine log grid and the 400x400 map is then
built two ways -- pixel CENTRE (reproducing the pipeline) and pixel AVERAGED
(supersampled, the correct cell average) -- and both are pushed through the
Stage 4 concentration statistic, the figure's 3 deg top-hat, and (against the
G1.1 view suite) the Pearson correlation.

Cheap: no 400x400xn_s array is ever formed, so it runs in seconds in a few
hundred MB. Needs numpy and scipy only.
"""
from __future__ import annotations

import sys
import numpy as np

from pathlib import Path                                          # noqa: E402
# The gNFW^2 line-of-sight integrator is vendored (gnfw_jmap.py, pure numpy),
# so this runs without the separate gce-fisher project on the path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gnfw_jmap import j_map, DEFAULTS as A2020_DEFAULTS           # noqa: E402
from scipy.ndimage import gaussian_filter                          # noqa: E402
from scipy.signal import fftconvolve                               # noqa: E402

NPIX, STEP = 400, 0.1
CTR = np.linspace(-20 + STEP / 2, 20 - STEP / 2, NPIX)
L2, B2 = np.meshgrid(CTR, CTR)
THETA = np.hypot(L2, B2)                    # stage4_stats.py aperture metric
MAPS = "/private/nfs/Data2/hestia_muru/stage3_maps"


def psi_of(l_deg, b_deg):
    """Great-circle angle from the GC -- what j_map's cos_psi uses."""
    l, b = np.deg2rad(l_deg), np.deg2rad(b_deg)
    return np.rad2deg(np.arccos(np.clip(np.cos(b) * np.cos(l), -1, 1)))


def radial_profile(n_psi=4000, n_s=3200):
    """J(psi) on a log grid, with the LOS convergence printed."""
    psi = np.exp(np.linspace(np.log(1e-5), np.log(35.0), n_psi))
    ref = j_map(psi, np.array([0.0]), n_s=800)[0]
    for n in (400, n_s):
        d = np.abs(j_map(psi, np.array([0.0]), n_s=n)[0] / ref - 1.0).max()
        print(f"  LOS n_s={n:5d} vs 800: max rel diff {d:.3e}")
    return psi, j_map(psi, np.array([0.0]), n_s=n_s)[0]


def build_maps():
    psi, J = radial_profile()
    lo, lJ = np.log(psi), np.log(J)

    def Jof(p):
        return np.exp(np.interp(np.log(np.clip(p, psi[0], psi[-1])), lo, lJ))

    def cell_average(nsub, lsel, bsel):
        off = (np.arange(nsub) + 0.5) / nsub - 0.5
        du, dv = np.meshgrid(off * STEP, off * STEP)
        acc = np.zeros(lsel.shape)
        for a, b in zip(du.ravel(), dv.ravel()):
            acc += Jof(psi_of(lsel + a, bsel + b))
        return acc / (nsub * nsub)

    centre = Jof(psi_of(L2, B2))
    # 9x9 everywhere, refined toward the cusp where the midpoint error lives
    avg = cell_average(9, L2, B2)
    for half, nsub in ((2.0, 61), (0.35, 401)):
        sel = (np.abs(L2) < half) & (np.abs(B2) < half)
        avg[sel] = cell_average(nsub, L2[sel], B2[sel])
    return centre, avg


def concentration(m, fwhm, r_in=3.0):
    s = gaussian_filter(m, fwhm / 2.3548 / STEP)
    return s[THETA < r_in].sum() / s[THETA < 10].sum()


def concentration_tophat(m, r_in=3.0):
    """The FIGURE's beam: 3 deg flat-(l,b) top-hat, figure_render.py:178."""
    kern = np.hypot(*np.meshgrid(*(np.arange(-30, 31) * STEP,) * 2)) <= 3.0
    s = fftconvolve(m, kern / kern.sum(), mode="same")
    return s[THETA < r_in].sum() / s[THETA < 10].sum()


def main():
    print("A2020_DEFAULTS:", A2020_DEFAULTS)
    centre, avg = build_maps()

    print("\nC = F(<3 deg)/F(<10 deg) for the spherical gNFW^2 template")
    print("  smoothing              pixel-centre   pixel-averaged")
    for label, fwhm in (("Gaussian FWHM 2 deg", 2.0), ("Gaussian FWHM 3 deg", 3.0)):
        print(f"  {label:22s} {concentration(centre, fwhm):.4f}         "
              f"{concentration(avg, fwhm):.4f}")
    print(f"  {'3 deg flat top-hat':22s} {concentration_tophat(centre):.4f}     "
          f"    {concentration_tophat(avg):.4f}")

    print(f"\n  central-pixel ratio averaged/centre: "
          f"{avg[199, 199] / centre[199, 199]:.3f}")
    print(f"  fraction of F(<3 deg) inside 0.1 deg (unsmoothed, averaged): "
          f"{avg[THETA < 0.1].sum() / avg[THETA < 3].sum():.3f}")

    # correlation shift, measured on ONE halo (36 views) -- see the note
    try:
        z = np.load(f"{MAPS}/G1.1.npz")
    except FileNotFoundError:
        print("\n  (stage3_maps not mounted; skipping the correlation check)")
        return
    print("\nPearson shift on G1.1's 36 views")
    for fwhm in (2.0, 3.0):
        sig = fwhm / 2.3548 / STEP
        tc = gaussian_filter(centre, sig).ravel()
        ta = gaussian_filter(avg, sig).ravel()
        for key in ("jmap", "starsold"):
            rc, ra = [], []
            for v in range(len(z["azimuth_deg"])):
                m = gaussian_filter(np.asarray(z[key][v], float), sig).ravel()
                rc.append(np.corrcoef(m, tc)[0, 1])
                ra.append(np.corrcoef(m, ta)[0, 1])
            rc, ra = np.array(rc), np.array(ra)
            print(f"  FWHM {fwhm}  {key:9s}  median {np.median(rc):.4f} -> "
                  f"{np.median(ra):.4f}   max |shift| {np.abs(ra - rc).max():.4f}")


if __name__ == "__main__":
    main()
