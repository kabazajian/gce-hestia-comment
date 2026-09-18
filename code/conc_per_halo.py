"""Per-halo concentration table + sensitivity to the unresolved centre.

Oscar Macias (2026-09-16): make explicit what C = 0.43 summarises, give the
per-halo range over the 36 azimuths, and show how sensitive the result is to
the treatment of the unresolved centre, since the numerator F(<3 deg) includes
the inner degrees the simulation does not resolve (softening 0.22 kpc = 1.6 deg
at the solar circle).

The sensitivity test: an ANNULAR concentration that excises the core from both
numerator and denominator,

    C(r_cut) = F(r_cut <= theta < 3 deg) / F(r_cut <= theta < 10 deg),

for r_cut = 0 (the published C), 1.0, 1.6 (the softening subtense) and 2.0 deg,
on the same FWHM-smoothed maps stage4_stats.py used. If the ranking
J > every bulge template > stars survives r_cut = 1.6-2.0 deg, the result does
not rest on the unresolved centre. Also reports the fraction of F(<3 deg) that
lies inside 1.6 deg, i.e. how much of the numerator is unresolved.

Templates: F98 / Cao13 / Coleman20 from stage4_templates.npz; gNFW^2 is the
CELL-AVERAGED template (check_gnfw2_quadrature.build_maps), the corrected
0.51 version, not the pixel-centre 0.49 one. Light: one halo (~90 MB) in
memory at a time.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

DATA = Path("/private/nfs/Data2/hestia_muru")
OUT = Path(__file__).resolve().parents[1] / "figures"
HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
NPIX, STEP = 400, 0.1
CTR = np.linspace(-20 + STEP / 2, 20 - STEP / 2, NPIX)
L2, B2 = np.meshgrid(CTR, CTR)
THETA = np.hypot(L2, B2)                     # stage4_stats.py aperture metric
FWHMS = (2.0, 3.0)
RCUTS = (0.0, 1.0, 1.6, 2.0)
TEMPLATES = ("f98", "cao13", "coleman20", "gnfw2")


def smooth(m, fwhm):
    return gaussian_filter(np.asarray(m, float), fwhm / 2.3548 / STEP)


def conc(m, r_cut, r_in=3.0, r_out=10.0):
    num = m[(THETA >= r_cut) & (THETA < r_in)].sum()
    den = m[(THETA >= r_cut) & (THETA < r_out)].sum()
    return float(num / den)


def core_fraction(m, r_core=1.6, r_in=3.0):
    """Fraction of the numerator F(<3 deg) that lies inside r_core."""
    return float(m[THETA < r_core].sum() / m[THETA < r_in].sum())


def load_templates():
    z = np.load(DATA / "stage4_templates.npz")
    t = {k: z[k].astype(float) for k in ("f98", "cao13", "coleman20")}
    from check_gnfw2_quadrature import build_maps
    _, t["gnfw2"] = build_maps()             # cell-averaged, the 0.51 version
    return t


def main():
    tm = load_templates()
    tsm = {f: {k: smooth(v, f) for k, v in tm.items()} for f in FWHMS}
    res = {"rcuts_deg": RCUTS, "fwhms_deg": FWHMS, "templates": {}, "halos": {}}
    for f in FWHMS:
        res["templates"][str(f)] = {
            k: {str(r): conc(tsm[f][k], r) for r in RCUTS} for k in TEMPLATES}

    allJ = {f: {r: [] for r in RCUTS} for f in FWHMS}
    allS = {f: {r: [] for r in RCUTS} for f in FWHMS}
    for h in HALOS:
        z = np.load(DATA / "stage3_maps" / f"{h}.npz")
        n = len(z["azimuth_deg"])
        rec = {str(f): {} for f in FWHMS}
        for f in FWHMS:
            J = [smooth(z["jmap"][v], f) for v in range(n)]
            S = [smooth(z["starsold"][v], f) for v in range(n)]
            for r in RCUTS:
                cj = np.array([conc(m, r) for m in J])
                cs = np.array([conc(m, r) for m in S])
                allJ[f][r] += cj.tolist(); allS[f][r] += cs.tolist()
                rec[str(f)][str(r)] = {
                    "J_median": float(np.median(cj)), "J_min": float(cj.min()),
                    "J_max": float(cj.max()),
                    "stars_median": float(np.median(cs)),
                    "stars_min": float(cs.min()), "stars_max": float(cs.max()),
                    "J_gt_own_stars": int((cj > cs).sum()),
                    "J_gt_template": {k: int((cj > res["templates"][str(f)][k][str(r)]).sum())
                                      for k in TEMPLATES},
                }
            rec[str(f)]["core_fraction_1p6_J"] = float(np.median([core_fraction(m) for m in J]))
            rec[str(f)]["core_fraction_1p6_stars"] = float(np.median([core_fraction(m) for m in S]))
            del J, S
        rec["n_views"] = n
        res["halos"][h] = rec
        del z
        print(f"{h}: done", flush=True)

    res["all216"] = {}
    for f in FWHMS:
        res["all216"][str(f)] = {}
        for r in RCUTS:
            cj, cs = np.array(allJ[f][r]), np.array(allS[f][r])
            res["all216"][str(f)][str(r)] = {
                "J_median": float(np.median(cj)), "J_min": float(cj.min()), "J_max": float(cj.max()),
                "stars_median": float(np.median(cs)), "stars_min": float(cs.min()), "stars_max": float(cs.max()),
                "J_gt_own_stars": int((cj > cs).sum()),
                "J_gt_template": {k: int((cj > res["templates"][str(f)][k][str(r)]).sum()) for k in TEMPLATES},
            }
    out = OUT / "conc_per_halo_2026-09-16.json"
    out.write_text(json.dumps(res, indent=1))

    # ---- print -----------------------------------------------------------
    for f in FWHMS:
        print(f"\n===== FWHM {f} deg =====")
        print(f"{'r_cut':>5s} | {'J med [min-max]':>22s} | {'stars med [min-max]':>22s} | "
              f"{'F98':>5s} {'Cao':>5s} {'Col':>5s} {'gNFW2':>5s} | J>F98 J>Cao J>Col J>gNFW2 J>stars")
        for r in RCUTS:
            a = res["all216"][str(f)][str(r)]; t = res["templates"][str(f)]
            g = a["J_gt_template"]
            print(f"{r:5.1f} | {a['J_median']:.3f} [{a['J_min']:.3f}-{a['J_max']:.3f}] | "
                  f"{a['stars_median']:.3f} [{a['stars_min']:.3f}-{a['stars_max']:.3f}] | "
                  f"{t['f98'][str(r)]:.3f} {t['cao13'][str(r)]:.3f} {t['coleman20'][str(r)]:.3f} {t['gnfw2'][str(r)]:.3f} | "
                  f"{g['f98']:5d} {g['cao13']:5d} {g['coleman20']:5d} {g['gnfw2']:7d} {a['J_gt_own_stars']:7d}")
        print(f"\n  per halo, r_cut = 0 (published C) and 1.6 deg:")
        print(f"  {'halo':5s} {'C_J med [min-max]':>22s} {'C_stars med':>11s} | {'C_J(1.6) med [min-max]':>25s} {'C_st(1.6)':>9s} | core frac J / stars")
        for h in HALOS:
            a0 = res["halos"][h][str(f)]["0.0"]; a1 = res["halos"][h][str(f)]["1.6"]
            cf = res["halos"][h][str(f)]
            print(f"  {h:5s} {a0['J_median']:.3f} [{a0['J_min']:.3f}-{a0['J_max']:.3f}] {a0['stars_median']:11.3f} | "
                  f"{a1['J_median']:.3f} [{a1['J_min']:.3f}-{a1['J_max']:.3f}] {a1['stars_median']:9.3f} | "
                  f"{cf['core_fraction_1p6_J']:.2f} / {cf['core_fraction_1p6_stars']:.2f}")
    print("\nwrote", out)


if __name__ == "__main__":
    main()
