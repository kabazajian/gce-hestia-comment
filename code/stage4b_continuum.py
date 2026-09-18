"""Stage 4b: continuum test + typicality (plan §4b; Kev's hold-out design).

Light local analysis on stage4_stats.json (fetched from nucosmo). Three
questions, all answered as COUNTING statements (Ryan's framing rule — no
sigma, no CI, no population of Milky Ways):

1. PREFERENCE COUNTS — per view: is the corrected J-map better correlated
   with spherical gNFW2 or with each bulge combo? Reported as "N of 216
   views (M of 6 halos)".

2. CONTINUUM (hold-out) — does the 216-view statistic-vector cloud form one
   describable distribution? Candidate: multivariate Gaussian fit to the
   training views; a held-out view "fits" if its Mahalanobis d2 is inside
   the chi2_k 90% quantile. Rungs of increasing strictness:
     a. random ~10% (x20 repeats)          — azimuthal smoothness only
     b. contiguous 40-deg arcs, 1/halo     — unseen directions
     c. leave-one-halo-out (6-fold)        — the load-bearing rung
     d. leave-one-realization-out (3-fold) — strictest (pairs share a LG)
   Pass/fail per rung: held-out coverage vs the ~90% expected, counted.

3. TYPICALITY OF THE TEMPLATES — each template map's statistic vector
   ranked in the empirical 216-view distribution, per statistic
   (percentile = counting), plus its Mahalanobis quantile under the fitted
   Gaussian IF rung c/d validated it (otherwise per-halo counting only).

Statistic vector (FWHM 2 deg): corr_gnfw2, corr_combo_1, q30, c4_30,
log10 conc_3_10. Five dims for 216 (correlated) views; chosen before
looking at the hold-out outcomes.

Mirror check: reports max |stat(mirror) - stat| per statistic over G3.2's
views; the effective view count claimed in the write-up follows from it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

STATS = Path(sys.argv[1] if len(sys.argv) > 1
             else "/private/nfs/Data2/hestia_muru/stage4_stats.json")
FWHM = sys.argv[2] if len(sys.argv) > 2 else "2.0"   # Stage 5 smoothing axis
HALOS = ("G1.1", "G1.2", "G2.1", "G2.2", "G3.1", "G3.2")
PAIRS = (("G1.1", "G1.2"), ("G2.1", "G2.2"), ("G3.1", "G3.2"))
VEC = ("corr_gnfw2", "corr_combo_1", "q30", "c4_30", "log_conc_3_10")
CHI2_90 = 9.236   # chi2, k=5, 0.90 quantile


def vec_of(st):
    return np.array([st["corr_gnfw2"], st["corr_combo_1"], st["q30"],
                     st["c4_30"], np.log10(st["conc_3_10"])])


def main():
    d = json.load(open(STATS))
    combos = [t for t in d["_provenance"]["templates"] if t.startswith("combo_")]

    X, halo_of = [], []
    pref = {c: 0 for c in combos}
    pref_halos = {c: set() for c in combos}
    for h in HALOS:
        for row in d[h]:
            st = row[FWHM]["j"]
            X.append(vec_of(st))
            halo_of.append(h)
            for c in combos:
                if st["corr_gnfw2"] > st[f"corr_{c}"]:
                    pref[c] += 1
                    pref_halos[c].add(h)
    X = np.array(X); halo_of = np.array(halo_of)
    n = len(X)

    print(f"== 1. preference counts ({n} views, FWHM {FWHM} deg) ==")
    for c in combos:
        print(f"  corr(J, gnfw2) > corr(J, {c}): {pref[c]} of {n} views "
              f"({len(pref_halos[c])} of 6 halos)")

    print("\n== mirror check (G3.2) ==")
    for k in ("corr_gnfw2", "corr_combo_1", "q30", "c4_30"):
        dmax = max(abs(m[FWHM][k] - r[FWHM]["j"][k])
                   for m, r in zip(d["mirror_check_G3.2"], d["G3.2"]))
        print(f"  max |mirror - original| {k}: {dmax:.4f}")

    print("\n== 2. continuum hold-out (Gaussian, 90% Mahalanobis region) ==")

    def fit(train):
        mu = train.mean(0)
        C = np.cov(train.T)
        Ci = np.linalg.inv(C)
        return mu, Ci

    def inside(test, mu, Ci):
        dv = test - mu
        return np.einsum("ij,jk,ik->i", dv, Ci, dv) <= CHI2_90

    rng = np.random.default_rng(0)
    hits = tot = 0
    for _ in range(20):
        idx = rng.permutation(n)
        te, tr = idx[:22], idx[22:]
        mu, Ci = fit(X[tr])
        hits += inside(X[te], mu, Ci).sum(); tot += len(te)
    print(f"  a. random 10% x20: {hits} of {tot} held-out inside 90% "
          f"(expect ~{0.9*tot:.0f})")

    hits = tot = 0
    for h in HALOS:
        w = np.where(halo_of == h)[0]
        a0 = rng.integers(0, 36)
        arc = w[[(a0 + k) % 36 for k in range(4)]]
        tr = np.setdiff1d(np.arange(n), arc)
        mu, Ci = fit(X[tr])
        hits += inside(X[arc], mu, Ci).sum(); tot += len(arc)
    print(f"  b. 40-deg arcs (1/halo): {hits} of {tot} inside 90% "
          f"(expect ~{0.9*tot:.0f})")

    for label, groups in (("c. leave-one-halo-out", [(h,) for h in HALOS]),
                          ("d. leave-one-realization-out", PAIRS)):
        hits = tot = 0
        per = []
        for g in groups:
            te = np.isin(halo_of, g)
            mu, Ci = fit(X[~te])
            k = inside(X[te], mu, Ci).sum()
            per.append(f"{'+'.join(g)}:{k}/{te.sum()}")
            hits += k; tot += te.sum()
        print(f"  {label}: {hits} of {tot} inside 90% (expect ~{0.9*tot:.0f})")
        print(f"     per fold: {'  '.join(per)}")

    print("\n== 3. template typicality (percentile rank in the 216 views, "
          "per statistic; counting) ==")
    mu, Ci = fit(X)
    hdr = "  {:<12s}" + "".join(f" {v:>14s}" for v in VEC) + "  Mahal-q"
    print(hdr.format("template"))
    from scipy.stats import chi2
    for t in ("gnfw2", *combos, "f98", "coleman20", "nb"):
        v = vec_of(d["templates"][t][FWHM])
        ranks = [(X[:, i] < v[i]).mean() for i in range(len(VEC))]
        dv = v - mu
        mq = chi2.cdf(dv @ Ci @ dv, df=len(VEC))
        print(("  {:<12s}" + "".join(f" {r:>14.3f}" for r in ranks)
               + f"  {mq:>7.4f}").format(t))

    print("\n== J vs own old stars, per view ==")
    cj = [row[FWHM]["corr_j_starsold"] for h in HALOS for row in d[h]]
    print(f"  corr(J, starsold): min {min(cj):.3f} median {np.median(cj):.3f} "
          f"max {max(cj):.3f}")
    # which template does each map PREFER — J vs its own stellar map
    both = 0
    for h in HALOS:
        for row in d[h]:
            j, s = row[FWHM]["j"], row[FWHM]["starsold"]
            if (j["corr_gnfw2"] > j["corr_combo_1"]) and \
               (s["corr_combo_1"] > s["corr_gnfw2"]):
                both += 1
    print(f"  views where J prefers gnfw2 AND its own stars prefer the bulge "
          f"combo: {both} of {n}")


if __name__ == "__main__":
    main()
