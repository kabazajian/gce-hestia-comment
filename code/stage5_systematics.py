"""Stage 5 of validation/muru_redo_plan_2026-08-30.md: the systematics budget.

Runs ON nucosmo. Recomputes the Stage 4 statistic vector (corr_gnfw2,
corr_combo_1, corr_f98, q30, c4_30, conc_3_10 at FWHM 2 deg) for the
J-maps of all 36 views of G3.2 (the anchored halo) under each variation,
against an in-script baseline computed through the identical code path:

  k16 / k64     — Stage 2 kNN density recomputed at k=16 / k=64
  smax50        — line-of-sight truncation 50 kpc instead of 100
  rsun8p0       — observer at Muru's 8.0 kpc instead of our 8.18
  half_a/half_b — random half-samples of the DM particles, rho recomputed
                  per half (particle noise, the honest way)

plus, on G3.1 (the worst Stage 1 center shift, 0.27 kpc):

  center_ahf    — AHF center used raw instead of shrinking-spheres refined

The smoothing-scale axis (FWHM 2 vs 3 deg) needs no new maps — both are in
stage4_stats.json — and is analyzed locally.

Output: stage5_systematics.json beside the data — per variant, per view,
the statistic vector; the local summarizer turns it into max/median |delta|
and preference-flip counts.
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

import stage3_views as s3            # noqa: E402
import stage4_stats as s4            # noqa: E402
from particle_projection import knn_density  # noqa: E402

FWHM = 2.0
TEMPLATES = ("gnfw2", "combo_1", "f98")
AZ = np.deg2rad(np.arange(0.0, 360.0, 10.0))


def jmap_for(pos_c, mass, rho, z_hat, obs, s_max):
    l, b, s = s3.lbs(pos_c, obs, z_hat)
    sel = (np.abs(l) < 20) & (np.abs(b) < 20) & (s <= s_max)
    return s3.binmap(l[sel], b[sel], (mass * rho)[sel] / s[sel] ** 2) / s3.PIX_SR


def stats_all_views(pos_c, mass, rho, z_hat, x_hat, tmpl_sm,
                    r_sun=8.18, s_max=100.0):
    y_hat = np.cross(z_hat, x_hat)
    rows = []
    for a in AZ:
        obs = r_sun * (np.cos(a) * x_hat + np.sin(a) * y_hat)
        m = jmap_for(pos_c, mass, rho, z_hat, obs, s_max)
        rows.append(s4.stats_for(m, tmpl_sm, FWHM))
    return rows


def main():
    tz = np.load(DATA / "stage4_templates.npz")
    tmpl_sm = {n: s4.smooth(tz[n], FWHM) for n in TEMPLATES}
    out = {"_provenance": {"script": "scripts/muru_redo/stage5_systematics.py",
                           "date": "2026-08-30", "fwhm": FWHM,
                           "templates": TEMPLATES}}
    t0 = time.time()

    # ---- G3.2 variants ------------------------------------------------------
    pos_c, mass, rho32, star, old, z_hat, x_hat, fr = s3.load_halo("G3.2")
    dm = rho32 > 0                      # rho set only on dm rows
    print(f"G3.2 loaded [{time.time()-t0:.0f} s]", flush=True)

    out["base"] = stats_all_views(pos_c, mass, rho32, z_hat, x_hat, tmpl_sm)
    print(f"base [{time.time()-t0:.0f} s]", flush=True)

    for k in (16, 64):
        rk = np.zeros_like(rho32)
        rk[dm] = knn_density(pos_c[dm], mass[dm], k=k)
        out[f"k{k}"] = stats_all_views(pos_c, mass, rk, z_hat, x_hat, tmpl_sm)
        print(f"k{k} [{time.time()-t0:.0f} s]", flush=True)

    out["smax50"] = stats_all_views(pos_c, mass, rho32, z_hat, x_hat, tmpl_sm,
                                    s_max=50.0)
    out["rsun8p0"] = stats_all_views(pos_c, mass, rho32, z_hat, x_hat, tmpl_sm,
                                     r_sun=8.0)
    print(f"smax50 + rsun8p0 [{time.time()-t0:.0f} s]", flush=True)

    rng = np.random.default_rng(0)
    half = rng.random(dm.sum()) < 0.5
    for tag, sel_half in (("half_a", half), ("half_b", ~half)):
        idx = np.where(dm)[0][sel_half]
        keep = np.zeros(len(mass), bool)
        keep[idx] = True
        rh = np.zeros_like(rho32)
        rh[keep] = knn_density(pos_c[keep], mass[keep], k=32)
        out[tag] = stats_all_views(pos_c[keep], mass[keep], rh[keep],
                                   z_hat, x_hat, tmpl_sm)
        print(f"{tag} [{time.time()-t0:.0f} s]", flush=True)

    # ---- G3.1 center variant ------------------------------------------------
    pos_c, mass, rho, star, old, z_hat, x_hat, fr = s3.load_halo("G3.1")
    out["g31_base"] = stats_all_views(pos_c, mass, rho, z_hat, x_hat, tmpl_sm)
    off = np.array(json.load(open(DATA / "muru_frames.json"))
                   ["G3.1"]["center_refined_offset_kpc"])
    out["g31_center_ahf"] = stats_all_views(pos_c + off, mass, rho,
                                            z_hat, x_hat, tmpl_sm)
    print(f"G3.1 center pair [{time.time()-t0:.0f} s]", flush=True)

    with open(DATA / "stage5_systematics.json", "w") as f:
        json.dump(out, f)
    print("wrote", DATA / "stage5_systematics.json", flush=True)


if __name__ == "__main__":
    main()
