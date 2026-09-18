"""Stage 4 template construction (validation/muru_redo_plan_2026-08-30.md §4).

Builds the comparison-template maps the Stage 4 statistics run on nucosmo
needs, and writes them to the data directory over NFS. Runs LOCALLY on the
Studio because it needs repo data (Gordon/Pohl FITS, gcepy, the ranking
artifacts); it is seconds of compute.

Templates, all pure morphology on the 0.1 deg analysis grid, gcepy
orientation [ib, il] (l, b increasing with index — simulated maps match):

  gnfw2       — spherical gNFW^2 gamma=1.2 J-map, dm_templates.j_map at
                A2020_DEFAULTS (R_sun = 8.18) — no PSF, no exposure.
  f98, coleman20 — Gordon repository FITS (the maps Coleman et al. 2020 cite),
                flux-conserving regrid 0.2 -> 0.1 deg, flipped to gcepy
                orientation. Orientation is CHECKED against gcepy's own
                energy-summed templates (corr correct > corr mirrored), not
                assumed — this project has met seven handedness traps.
  nb          — Pohl 2022 Nuclear_Bulge (1 arcmin, offset CRVAL — regrid_car
                reads the WCS, which handles it).
  combo_1..N  — the "lambda-superimposed" bulge models (Kev, 2026-08-30):
                for each of the TOP_N (bulge+NB, diffuse-model) cells of the
                main-analysis rankings by ln_H, the map
                w_bulge * bulge_hat + w_nb * nb_hat, with w = the cell's
                FITTED photon counts summed over energy bins
                (sum_b 10^theta[b,c] * template_total[b,c]) and hat = unit-sum
                morphology. Cells whose bulge is cao13 use gcepy's
                energy-summed cao13 (PSF+exposure baked — stated caveat;
                no pure-morphology cao13 map exists in this repo).

Caveat recorded in the output metadata: fitted counts weight the baked
(PSF x exposure) templates, and we apply those weights to pure morphologies —
exact for the ratio's meaning (photons per component), approximate in shape
by the smooth exposure gradient; irrelevant at the 2-3 deg smoothing of
every downstream statistic.

The A42 bulge-lambda composite is the eventual "lambda superimposed" referent
and is NOT yet in the rankings (results/gnfw2_vs_bulge_ranking_summary_
2026-08-23.md, final caveat); these combos stand in, per the plan.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

OUT = Path("/private/nfs/Data2/hestia_muru")
TOP_N = 10
NPIX, STEP = 400, 0.1
CTR = np.linspace(-20 + STEP / 2, 20 - STEP / 2, NPIX)   # pixel centers

RANK_FILES = ["runs/gnfw2_vs_bulge_evidence.json",
              "runs/gnfw2_vs_bulge_evidence_galp21.json"]


#: 🚨 ln_H is comparable only WITHIN a data selection: cells that mask more
#: pixels have fewer data and trivially higher ln_H, so a pooled ranking is a
#: ranking of data volume, not of models (the first draft of this script did
#: exactly that and returned ten rfree0_mask90 rows). The combos are therefore
#: ranked within ONE complete-data cell: rfree10_nomask, no plane mask — the
#: fullest dataset, inner point sources free (and, per the ranking summary,
#: the cell where the bulge-vs-gNFW2 comparison is cleanest). Recorded as a
#: choice, revisitable.
RANK_CELL = ("rfree10_nomask", None)


def top_cells():
    rows = []
    for f in RANK_FILES:
        rows += json.load(open(ROOT / f))["results"]
    ok = [r for r in rows
          if r["ok"] and r["hypothesis"].endswith("+NB")
          and r["ps_masking"] == RANK_CELL[0] and r["plane_mask"] == RANK_CELL[1]
          and r["ln_H"] is not None]
    ok.sort(key=lambda r: r["ln_H"], reverse=True)
    return ok[:TOP_N]


def fitted_weights(row, totals):
    """Energy-summed fitted counts of (bulge, nb) — the LAST TWO theta columns
    (names = diffuse + (bubble, isotropic, ps_fixed, ext_fixed) + excess, see
    gnfw2_vs_bulge_scan.py)."""
    bulge = row["hypothesis"].split("+")[0]
    th = np.array(row["theta"])                      # [n_bins, n_comp], log10
    wb = float(np.sum(10.0 ** th[:, -2] * totals[bulge]))
    wn = float(np.sum(10.0 ** th[:, -1] * totals["nb"]))
    return bulge, wb, wn


def load_morph_maps():
    """Pure-morphology maps in gcepy orientation, plus gcepy references."""
    from astropy.io import fits
    from gcereduce.regrid import regrid_car
    from gcefisher import likelihood as L

    maps, checks = {}, {}
    # analysis grid in FITS/WCS convention (matches the gtbin CCUBE:
    # CONVENTIONS.md § gtbin CCUBE longitude axis) — l DEcreasing with column;
    # the [:, ::-1] below flips to gcepy orientation.
    tgt = fits.Header()
    tgt["NAXIS"] = 2
    tgt["NAXIS1"] = tgt["NAXIS2"] = NPIX
    tgt["CTYPE1"], tgt["CTYPE2"] = "GLON-CAR", "GLAT-CAR"
    tgt["CRVAL1"] = tgt["CRVAL2"] = 0.0
    tgt["CRPIX1"] = tgt["CRPIX2"] = NPIX / 2 + 0.5
    tgt["CDELT1"], tgt["CDELT2"] = -STEP, STEP

    for name, path in (
            ("f98", "data/external/Gordon_bulge_templates/F98_BoxyBulge_arxv1611.06644_Normalized.fits"),
            ("coleman20", "data/external/Gordon_bulge_templates/Bulge_modulated_Coleman_etal_2019_Normalized.fits"),
            ("nb", "data/external/Pohl_etal_2022/Nuclear_Bulge_Normalized.fits")):
        with fits.open(ROOT / path) as h:
            m = regrid_car(h[0].data.astype(float), h[0].header, tgt)
        maps[name] = m[:, ::-1]                       # FITS/WCS -> gcepy l order
        # orientation check against gcepy's own (energy-summed) template where
        # one exists; NB has none in gcepy — checked via f98/coleman20 instead
        if name != "nb":
            ref = np.asarray(L._template("low", name)).reshape(14, NPIX, NPIX).sum(0)
            r_ok = np.corrcoef(maps[name].ravel(), ref.ravel())[0, 1]
            r_mir = np.corrcoef(maps[name][:, ::-1].ravel(), ref.ravel())[0, 1]
            checks[name] = {"corr_correct": float(r_ok), "corr_mirrored": float(r_mir)}
            assert r_ok > r_mir, f"{name}: orientation check FAILED ({r_ok} vs {r_mir})"

    from gcereduce.dm_templates import j_map
    maps["gnfw2"] = j_map(CTR, CTR)                   # spherical — orientation-proof
    assert maps["gnfw2"].shape == (NPIX, NPIX)

    cao = np.asarray(L._template("low", "cao13")).reshape(14, NPIX, NPIX).sum(0)
    maps["cao13"] = cao                               # gcepy orientation already
    checks["cao13"] = "gcepy energy-summed counts (PSF+exposure baked) — no pure-morphology map exists"

    totals = {n: np.asarray(L._template("low", n)).reshape(14, -1).sum(1)
              for n in ("f98", "cao13", "coleman20")}
    totals["nb"] = np.load(ROOT / "runs/acb4f1769f61/nb_template.npy").reshape(14, -1).sum(1)
    return maps, checks, totals


def main():
    maps, checks, totals = load_morph_maps()
    unit = {k: v / v.sum() for k, v in maps.items()}

    cells = top_cells()
    meta = {"script": "scripts/muru_redo/stage4_templates.py",
            "date": "2026-08-30", "orientation_checks": checks,
            "note": ("combo weights = fitted photon counts (energy-summed) "
                     "applied to unit-sum pure morphologies; A42 lambda-"
                     "composite not yet ranked, these stand in"),
            "cells": []}
    out = {k: maps[k].astype(np.float32) for k in ("gnfw2", "f98", "coleman20",
                                                   "cao13", "nb")}
    for i, row in enumerate(cells, 1):
        bulge, wb, wn = fitted_weights(row, totals)
        combo = wb * unit[bulge] + wn * unit["nb"]
        out[f"combo_{i}"] = (combo / combo.sum()).astype(np.float32)
        meta["cells"].append({
            "rank": i, "model": row["model"], "family": row["family"],
            "hypothesis": row["hypothesis"], "ps_masking": row["ps_masking"],
            "plane_mask": row["plane_mask"], "ln_H": row["ln_H"],
            "w_bulge_counts": wb, "w_nb_counts": wn,
            "nb_fraction": wn / (wb + wn),
            "cao13_caveat": bulge == "cao13"})
        print(f"combo_{i}: {row['hypothesis']:14s} {row['model']:>10s} "
              f"{row['ps_masking']:>22s} lnH {row['ln_H']:.0f}  "
              f"NB fraction {wn/(wb+wn):.3f}")

    np.savez_compressed(OUT / "stage4_templates.npz", **out)
    with open(OUT / "stage4_templates_meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print("orientation checks:", checks)
    print("wrote", OUT / "stage4_templates.npz")


if __name__ == "__main__":
    main()
