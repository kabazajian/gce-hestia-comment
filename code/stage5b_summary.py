"""Summary of the six-halo kNN k test (stage5b_k_allhalos.py), no particles.

Reads the stage-5b JSON and prints, at Gaussian FWHM 2 and 3 deg, the largest
change over all 216 views of

  C    = F(<3)/F(<10)
  C'   = F(1.6<th<3)/F(1.6<th<10)   (inner 1.6 deg, the softening, excluded)
  core = F(<1.6)/F(<3)

relative to k = 32 and over the full span k = 16 -> 64, and checks at each k
that every annihilation map is more concentrated than the old-star map of the
same view and than each bulge template (conc_per_halo.json).

The Comment quotes the FWHM 2 deg numbers: 0.007 in C over the full span
(0.004 relative to k = 32) and 0.004 in C'.

Keys in the JSON are particle-FILE labels; nothing here is quoted per halo.
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_DERIVED = HERE.parent / "derived"
if (_DERIVED / "stage5b_k_allhalos.json").exists():           # public release
    K_JSON = _DERIVED / "stage5b_k_allhalos.json"
    T_JSON = _DERIVED / "conc_per_halo.json"
else:
    K_JSON = HERE.parent / "figures" / "stage5b_k_allhalos_2026-11-17.json"
    T_JSON = HERE.parent / "figures" / "conc_per_halo_2026-09-16.json"

d = json.load(open(K_JSON))
tmpl = json.load(open(T_JSON))["templates"]
halos = [h for h in d if not h.startswith("_")]
BULGE = ("f98", "cao13", "coleman20")
CUT = {"C": "0.0", "Cp": "1.6"}


def arr(key, f, q):
    return np.array([v[f][q] for h in halos for v in d[h][key]])


for f in ("2.0", "3.0"):
    print(f"Gaussian FWHM {f} deg, {len(arr('k32', f, 'C'))} views")
    for q in ("C", "Cp", "core"):
        a16, a32, a64 = (arr(f"k{k}", f, q) for k in (16, 32, 64))
        print(f"  {q:<4s} max|k16-k32| {np.abs(a16 - a32).max():.4f}"
              f"  max|k64-k32| {np.abs(a64 - a32).max():.4f}"
              f"  max|k16-k64| {np.abs(a16 - a64).max():.4f}"
              f"  median k16-k64 {np.median(a16 - a64):+.4f}")
    for q in ("C", "Cp"):
        st = arr("stars_old", f, q)
        tmax = max(tmpl[f][t][CUT[q]] for t in BULGE)
        for k in (16, 32, 64):
            a = arr(f"k{k}", f, q)
            print(f"  {q:<4s} k{k}: annihilation min {a.min():.3f} median {np.median(a):.3f};"
                  f" above own stars in {int((a > st).sum())}/{a.size},"
                  f" above every bulge template (max {tmax:.3f}) in {int((a > tmax).sum())}/{a.size}")
