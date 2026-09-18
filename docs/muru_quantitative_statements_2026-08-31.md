# Which quantitative statements the Muru Comment can safely make

**2026-08-31.** Kev asked for a quantitative statement in the Comment —
specifically the correlation of the corrected J maps with the observed
stellar-bulge templates. Measured first (from `stage4_stats.json` plus a
new shape-only statistic), because **two of the three candidate
correlation statements would support their claim, not ours.**

## 🚨 Do NOT quote map correlations on the J side

**Raw Pearson on the smoothed maps** (median over 216 views, FWHM 2°/3°):

| | gNFW² | F98 | Cao13 | Coleman20 |
|---|---|---|---|---|
| corrected J | 0.920 / 0.963 | 0.906 / 0.934 | 0.916 / 0.933 | **0.930** / 0.949 |
| old stars | 0.723 / 0.796 | 0.944 / 0.948 | 0.977 / 0.978 | 0.960 / 0.961 |

The corrected J maps correlate with **everything** at r ≈ 0.90–0.96 —
Coleman20 (0.930) marginally *above* gNFW² (0.920) at FWHM 2°. All these
maps are centrally peaked, so linear-map Pearson is dominated by that
common structure (the Stage 5 finding, now quantified against the bulge
templates). Quoting r(J, Coleman20) = 0.93 in the Comment would hand the
Reply a sentence.

**Shape-only correlation** — each map divided by its own azimuthally
averaged radial profile over 0.5–15°, so only angular structure remains
(new here; both smoothings agree):

| | gNFW² | F98 | Cao13 | Coleman20 |
|---|---|---|---|---|
| corrected J | +0.11 | +0.65 | **+0.92** | +0.69 |
| old stars | +0.04 | +0.64 | +0.95 | +0.71 |

🚨 **The corrected J maps' angular structure correlates with Cao13's at
0.92 — essentially as well as the old stars do (0.95).** This is real and
must not be hidden: both the DM and the stellar bulge are flattened
*toward the disk plane in the same orientation*, so their angular patterns
match. (The gNFW² column is degenerate, not informative: a spherical
template has no angular structure to correlate with.)

**Consequence:** the corrected maps are *not* "morphologically unlike" the
bulge in the pattern sense. The defensible claim is about **degree and
radial profile**, which is also exactly what the July analysis predicted
would matter (`muru_projection_analysis.md` §5: "Even if the isophote
shapes were similar, the radial profiles are not, and a bin-by-bin
template fit keys on exactly that").

## ✅ The statement that IS robust: concentration

C ≡ F(<3°)/F(<10°), FWHM 2° (3° in parentheses):

| | C | J more concentrated |
|---|---|---|
| corrected J, 216 views | **0.427** median, range 0.34–0.50 (0.390) | — |
| F98 | 0.267 (0.254) | **216/216** |
| Cao13 | 0.269 (0.255) | **216/216** |
| Coleman20 | 0.292 (0.275) | **216/216** |
| gNFW² (γ=1.2) | 0.491 (0.453) | 6/216 (0/216) |
| old stars, 216 views | 0.231 (0.219) | — |

**Every one of the 216 corrected maps is more centrally concentrated than
every bulge template, at both smoothings**, and sits near the cuspy
spherical profile; the same simulations' old stellar populations land at
0.23, consistent with the bulge templates. This is the ρ²-weighting
signature stated quantitatively, and it is the discriminator a
template fit actually keys on.

Supporting, also robust (Stage 4/5): axis ratio q30 median 0.76 for the
corrected maps against 0.56 (Cao13), 0.62 (Coleman20), 0.72 (F98) —
rounder in 203/216 and 198/216 views respectively; and boxiness c4 above
F98's −0.024 in 216/216 views. ⚠ c4 does **not** discriminate against
Cao13/Coleman20, whose c4 is positive (+0.036, +0.026).

## Recommended framing for the manuscript

Lead with concentration; state the flattening as a difference of degree;
**preempt the pattern-correlation point explicitly** rather than letting a
Reply raise it — the corrected maps are flattened in the same sense as the
bulge, and saying so costs nothing while the concentration and
axis-ratio numbers carry the argument.
