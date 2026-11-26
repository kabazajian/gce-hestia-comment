"""Concentration of the Letter's own constructions versus the corrected maps,
on the footing of the SM figure (2026-11-24).

Question: how much does the projection correction itself
change the concentration contrast between annihilation and stars?

Maps: the adopted figure's inputs (figure_maps.npz, bestaz frames): ONE view
per galaxy, observer at 8 kpc, 15 kpc line of sight, 3-deg-radius TOP-HAT
aperture. This is NOT the 216-view, Gaussian FWHM 2 deg suite of the
concentration table; do not mix the two sets of numbers.

C = F(<3 deg)/F(<10 deg), pixels summed in the flat (l, b) plane.
Columns: linear DM [int rho s^2 ds], their stars (same statistic),
quadratic DM (linear DM squared), corrected stars [int rho_* ds],
corrected annihilation [int rho^2 ds].

Keys in the npz are particle-FILE labels; output is keyed by PAPER label.
"""
import json
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

from figure_render import PAPER_LABEL, X

# Public release layout: read and write ../derived. Otherwise the data volume.
_DERIVED = Path(__file__).resolve().parent.parent / "derived"
if (_DERIVED / "figure_maps.npz").exists():
    NPZ = _DERIVED / "figure_maps.npz"
    OUT = _DERIVED / "stage6_conc_figure_footing.json"
else:
    NPZ = Path("/private/nfs/Data2/hestia_muru/figure_maps.npz")
    OUT = Path(__file__).parent.parent / "figures" / "stage6_conc_figure_footing_2026-11-24.json"

z = np.load(NPZ)
L, B = np.meshgrid(X, X)
R = np.hypot(L, B)
kern = np.hypot(*np.meshgrid(*(np.arange(-30, 31) * 0.1,) * 2)) <= 3.0
kern = kern / kern.sum()


def conc(m):
    return float(m[R < 3].sum() / m[R < 10].sum())


cols = ("lin_dm", "their_stars", "quad_dm", "corr_stars", "j_corr")
out = {"_provenance": {"script": "analysis/code/stage6_conc_figure_footing.py",
                       "maps": "figure_maps.npz", "kernel": "3 deg top-hat", "views": "one per galaxy",
                       "los_kpc": 15, "labels": "paper"}}
for h in sorted(PAPER_LABEL, key=PAPER_LABEL.get):
    mu = z[f"{h}_mu_dm"].astype(float)
    st = z[f"{h}_mu_stars"].astype(float)
    j = fftconvolve(z[f"{h}_jmap_muru_geom"], kern, mode="same")
    sc = fftconvolve(z[f"{h}_stars_col"], kern, mode="same")
    v = [conc(m) for m in (mu, st, mu ** 2, sc, j)]
    out[PAPER_LABEL[h]] = dict(zip(cols, v))
    out[PAPER_LABEL[h]]["quad_over_their_stars"] = v[2] / v[1]
    out[PAPER_LABEL[h]]["j_over_corr_stars"] = v[4] / v[3]
    print(PAPER_LABEL[h], " ".join(f"{x:.3f}" for x in v),
          f"| quad/stars {v[2]/v[1]:.2f}  J/stars {v[4]/v[3]:.2f}")
a = np.array([[out[g][c] for c in cols] for g in out if g[0] == "G"])
out["median"] = dict(zip(cols, np.median(a, 0).tolist()))
print("median", np.round(np.median(a, 0), 3))
OUT.write_text(json.dumps(out, indent=1))
