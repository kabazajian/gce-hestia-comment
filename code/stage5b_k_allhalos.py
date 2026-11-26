"""Stage 5b: kNN k-sensitivity of the concentration, all six halos.

Runs ON nucosmo. Stage 5 varied k on one halo only (file label G3.2); this
stage tests the sensitivity at the centre on all six. For every halo and
each of k = 16, 32, 64, recompute the J-map of all 36 views and record, at
Gaussian FWHM 2 and 3 deg:

  C    = F(<3)/F(<10)                  (conc_per_halo.conc, r_cut = 0)
  C'   = F(1.6<th<3)/F(1.6<th<10)      (r_cut = 1.6 deg, the softening)
  core = F(<1.6)/F(<3)                 (conc_per_halo.core_fraction)

k = 32 is recomputed through the same code path rather than read from
rho_k32/, so the baseline differs from the variants in k alone. The old-star
maps do not depend on k and are recorded once per view.

Keys are particle-FILE labels (CLAUDE.md rule 6).

Output: stage5b_k_allhalos.json beside the data.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
DATA = Path(os.environ.get("MURU_DATA", "/Volumes/Data2/hestia_muru"))
sys.path.insert(0, str(DATA))

import stage3_views as s3            # noqa: E402
from particle_projection import knn_density  # noqa: E402

HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
KS = (16, 32, 64)
FWHMS = (2.0, 3.0)
AZ = np.deg2rad(np.arange(0.0, 360.0, 10.0))
NPIX, STEP = 400, 0.1
CTR = np.linspace(-20 + STEP / 2, 20 - STEP / 2, NPIX)
THETA = np.hypot(*np.meshgrid(CTR, CTR))


def smooth(m, fwhm):
    return gaussian_filter(np.asarray(m, float), fwhm / 2.3548 / STEP)


def conc(m, r_cut, r_in=3.0, r_out=10.0):
    num = m[(THETA >= r_cut) & (THETA < r_in)].sum()
    den = m[(THETA >= r_cut) & (THETA < r_out)].sum()
    return float(num / den)


def core_fraction(m, r_core=1.6, r_in=3.0):
    return float(m[THETA < r_core].sum() / m[THETA < r_in].sum())


def stats(m):
    out = {}
    for f in FWHMS:
        ms = smooth(m, f)
        out[str(f)] = {"C": conc(ms, 0.0), "Cp": conc(ms, 1.6),
                       "core": core_fraction(ms)}
    return out


def views(pos_c, w, z_hat, x_hat):
    """Stats of the 1/s^2-weighted map of weights w at each azimuth."""
    y_hat = np.cross(z_hat, x_hat)
    rows = []
    for a in AZ:
        obs = s3.R_SUN * (np.cos(a) * x_hat + np.sin(a) * y_hat)
        l, b, s = s3.lbs(pos_c, obs, z_hat)
        sel = (np.abs(l) < 20) & (np.abs(b) < 20) & (s <= s3.S_MAX)
        rows.append(stats(s3.binmap(l[sel], b[sel], w[sel] / s[sel] ** 2)
                          / s3.PIX_SR))
    return rows


def main():
    out = {"_provenance": {"script": "analysis/code/stage5b_k_allhalos.py",
                           "date": time.strftime("%Y-%m-%d"),
                           "ks": KS, "fwhms": FWHMS, "labels": "particle-file"}}
    t0 = time.time()
    for h in HALOS:
        pos_c, mass, rho32, star, old, z_hat, x_hat, fr = s3.load_halo(h)
        dm = rho32 > 0                      # rho set only on dm rows
        rec = {"stars_old": views(pos_c[old], mass[old], z_hat, x_hat)}
        for k in KS:
            rk = knn_density(pos_c[dm], mass[dm], k=k)
            if k == 32:
                rec["k32_vs_file_maxrel"] = float(
                    np.max(np.abs(rk / rho32[dm] - 1.0)))
            rec[f"k{k}"] = views(pos_c[dm], mass[dm] * rk, z_hat, x_hat)
            print(f"{h} k{k} [{time.time()-t0:.0f} s]", flush=True)
        out[h] = rec
        with open(DATA / "stage5b_k_allhalos.json", "w") as f:
            json.dump(out, f)
    print("wrote", DATA / "stage5b_k_allhalos.json", flush=True)


if __name__ == "__main__":
    main()
