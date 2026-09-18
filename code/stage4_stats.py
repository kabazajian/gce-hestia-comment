"""Stage 4 statistics pass (validation/muru_redo_plan_2026-08-30.md §4).

Runs ON nucosmo. For every Stage 3 view (216) and every Stage 4 template,
at Gaussian smoothing FWHM 2 and 3 deg (>= the 1.58 deg softening subtense):

  corr_<t>   — Pearson correlation with each template map (linear pixels,
               the analog of a linear template fit).
  q50/q30/q10 — isophote axis ratio at 50/30/10% of max: sqrt of the
               eigenvalue ratio of the second-moment matrix of the pixels
               above the level (about the GC), i.e. the moment analog of
               Muru et al.'s per-contour SVD.
  c4_50/c4_30 — boxiness: Fourier a4/a0 of the contour radius r(phi) in the
               major-axis frame. c4 < 0 = boxy, > 0 = disky (standard
               isophote convention).
  conc_1_10 / conc_3_10 — f(<1 deg)/f(<10 deg), f(<3)/f(<10).

Per view additionally:
  corr_j_starsold — J vs the SAME view's old-star map (their "DM and MSP
               morphologies are indistinguishable" claim, done with the
               correct functional on the DM side).
  stats for the old-star map itself (so the stellar side of the comparison
               goes through the identical statistic code).

Mirror check (G3.2, all views): every statistic recomputed on the
l-mirrored map. corr vs a TILTED template (F98 bar) is NOT mirror-invariant
a priori — measured, not assumed, and the result decides the effective
view count in Stage 4b.

Output: stage4_stats.json beside the data.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

DATA = Path(os.environ.get("MURU_DATA", "/Volumes/Data2/hestia_muru"))
MAPS = DATA / "stage3_maps"
FWHMS = (2.0, 3.0)
NPIX, STEP = 400, 0.1
CTR = np.linspace(-20 + STEP / 2, 20 - STEP / 2, NPIX)
L2, B2 = np.meshgrid(CTR, CTR)          # [ib, il]
THETA = np.hypot(L2, B2)
HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
NPHI, R_SCAN = 256, np.arange(0.2, 19.9, 0.05)


def smooth(m, fwhm):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(np.asarray(m, float), fwhm / 2.3548 / STEP)


def _contour_radius(m, level):
    """r(phi): outermost radius where the map still exceeds level*max."""
    from scipy.ndimage import map_coordinates
    phi = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)
    # rays in pixel coordinates about the GC pixel boundary (199.5, 199.5)
    il = 199.5 + np.outer(np.cos(phi), R_SCAN) / STEP
    ib = 199.5 + np.outer(np.sin(phi), R_SCAN) / STEP
    v = map_coordinates(m, [ib, il], order=1, mode="nearest")
    above = v >= level * m.max()
    idx = np.where(above.any(axis=1),
                   above.shape[1] - 1 - above[:, ::-1].argmax(axis=1), 0)
    return phi, R_SCAN[idx]


def isophote_stats(m, levels=(0.5, 0.3, 0.1), c4_levels=(0.5, 0.3)):
    out = {}
    for lv in levels:
        sel = m >= lv * m.max()
        x, y = L2[sel], B2[sel]
        C = np.cov(np.vstack([x, y]))
        ev, evec = np.linalg.eigh(C)
        out[f"q{int(lv*100)}"] = float(np.sqrt(max(ev[0], 0) / max(ev[1], 1e-30)))
        if lv in c4_levels:
            pa = np.arctan2(evec[1, 1], evec[0, 1])       # major axis
            phi, r = _contour_radius(m, lv)
            ph = phi - pa
            a0 = r.mean()
            c4 = 2.0 * np.mean(r * np.cos(4 * ph)) / a0
            out[f"c4_{int(lv*100)}"] = float(c4)
    return out


def stats_for(m_raw, tmpl_sm, fwhm):
    m = smooth(m_raw, fwhm)
    st = {f"corr_{n}": float(np.corrcoef(m.ravel(), t.ravel())[0, 1])
          for n, t in tmpl_sm.items()}
    st.update(isophote_stats(m))
    tot = m.sum()
    st["conc_1_10"] = float(m[THETA < 1].sum() / m[THETA < 10].sum())
    st["conc_3_10"] = float(m[THETA < 3].sum() / m[THETA < 10].sum())
    st["total"] = float(tot)
    return st


def main():
    t0 = time.time()
    tz = np.load(DATA / "stage4_templates.npz")
    tnames = list(tz.files)
    tmpl_sm = {f: {n: smooth(tz[n], f) for n in tnames} for f in FWHMS}

    out = {"_provenance": {
        "script": "scripts/muru_redo/stage4_stats.py", "date": "2026-08-30",
        "fwhms_deg": FWHMS, "templates": tnames,
        "c4_convention": "contour-radius a4/a0 in major-axis frame; <0 boxy"}}

    # templates through the identical code (correlations vs the OTHER
    # templates included — gives e.g. corr(combo_i, gnfw2) for free)
    out["templates"] = {n: {str(f): stats_for(tz[n], tmpl_sm[f], f)
                            for f in FWHMS} for n in tnames}
    print(f"templates done [{time.time()-t0:.0f} s]", flush=True)

    for halo in HALOS:
        z = np.load(MAPS / f"{halo}.npz")
        hv = []
        for i in range(z["jmap"].shape[0]):
            row = {"azimuth_deg": float(z["azimuth_deg"][i])}
            for f in FWHMS:
                j = smooth(z["jmap"][i], f)
                so = smooth(z["starsold"][i], f)
                row[str(f)] = {
                    "j": stats_for(z["jmap"][i], tmpl_sm[f], f),
                    "starsold": stats_for(z["starsold"][i], tmpl_sm[f], f),
                    "corr_j_starsold": float(
                        np.corrcoef(j.ravel(), so.ravel())[0, 1]),
                }
            hv.append(row)
        out[halo] = hv
        print(f"{halo}: {len(hv)} views [{time.time()-t0:.0f} s]", flush=True)

    # mirror check on G3.2: identical stats on the l-mirrored maps
    z = np.load(MAPS / "G3.2.npz")
    mir = []
    for i in range(z["jmap"].shape[0]):
        m = z["jmap"][i][:, ::-1]
        mir.append({str(f): stats_for(m, tmpl_sm[f], f) for f in FWHMS})
    out["mirror_check_G3.2"] = mir
    print(f"mirror check done [{time.time()-t0:.0f} s]", flush=True)

    with open(DATA / "stage4_stats.json", "w") as fh:
        json.dump(out, fh)
    print("wrote", DATA / "stage4_stats.json", flush=True)


if __name__ == "__main__":
    main()
