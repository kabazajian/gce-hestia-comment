# Per-halo concentration table and the unresolved-centre test (Oscar's ask)

⚠ **Relabelled 2026-09-16.** Per-halo rows now carry the PAPER's labels; the
particle files are transposed within pairs relative to Ref. [1] (Finding 7 in
`muru_ahf_frames_2026-09-15.md`). The all-216 rows, the templates, and every
ranking statement are unaffected. `conc_per_halo_2026-09-16.json` is still
keyed by particle-file label.

**2026-09-16.** Oscar Macias asked for a compact SM table saying what
C = 0.43 summarises, giving each halo's range over its 36 azimuths, making
"216 views of six halos" explicit, and — the substantive part — showing how
sensitive the result is to the unresolved centre, since the numerator
F(<3°) includes the inner degrees the simulation cannot resolve.

Script: `analysis/code/conc_per_halo.py` (light; one halo in memory at a
time). Output: `analysis/figures/conc_per_halo_2026-09-16.json`. Same maps,
same smoothing, same aperture metric as `stage4_stats.py`; the gNFW² entry
is the **cell-averaged** (0.51) template.

## The answer: the result does not rest on the unresolved centre

The test is an annular concentration that excises the core from both
numerator and denominator,
C(r_cut) = F(r_cut ≤ θ < 3°) / F(r_cut ≤ θ < 10°), for r_cut = 1°, **1.6°**
(the softening subtense), 2°.

| FWHM 2° | J median [range] | stars median | F98 | Cao | Coleman | gNFW² | J > every template | J > own stars |
|---|---|---|---|---|---|---|---|---|
| r_cut = 0 (published C) | 0.427 [0.342–0.496] | 0.231 | 0.267 | 0.269 | 0.292 | 0.508 | **216/216** | **216/216** |
| r_cut = 1.0° | 0.380 [0.304–0.441] | 0.205 | 0.241 | 0.240 | 0.262 | 0.417 | **216/216** | **216/216** |
| r_cut = 1.6° | 0.305 [0.243–0.353] | 0.165 | 0.198 | 0.194 | 0.213 | 0.300 | **216/216** | **216/216** |
| r_cut = 2.0° | 0.236 [0.188–0.273] | 0.128 | 0.156 | 0.151 | 0.168 | 0.215 | **216/216** | **216/216** |

Same at FWHM 3° (216/216 everywhere; raw numbers in the JSON). Every
concentration drops when the core is removed — but **no ranking changes at
any cut, at either smoothing.** The least concentrated view of the least
concentrated halo (G2.2, 0.243 at r_cut = 1.6°) still clears the most
concentrated bulge template (Coleman, 0.213).

**How much of the numerator is nominally unresolved:** 38–44% of F(<3°)
for J lies inside 1.6° (per-halo medians), 33–36% for the stars. State it
plainly; it is an *upper bound* on the unresolved contribution, since those
lines of sight also traverse well-resolved material at larger 3D radii.

⚠ **One thing that IS smoothing-dependent — do not oversell it.** With the
core excised, the simulated maps are as concentrated as the spherical
gNFW² cusp at FWHM 2° (J exceeds it in 118/216 views at 1.6°, 126/216 at
2°) but not at FWHM 3° (9/216, 20/216). Once the central spike is removed
the cusp loses its main advantage, so the near-equality is real, but it
sits within the smoothing systematic. Mention as a remark or leave out; it
is not load-bearing.

## Per halo (FWHM 2°): azimuthal spread is smaller than halo-to-halo spread

| halo | C_J median [range over 36 az.] | C_stars median | C_J(>1.6°) median [range] | C_stars(>1.6°) | core frac. J / stars |
|---|---|---|---|---|---|
| G1.1 | 0.456 [0.444–0.469] | 0.223 | 0.329 [0.322–0.339] | 0.160 | 0.41 / 0.33 |
| G1.2 | 0.384 [0.356–0.435] | 0.223 | 0.274 [0.252–0.311] | 0.159 | 0.39 / 0.34 |
| G2.1 | 0.350 [0.342–0.362] | 0.210 | 0.250 [0.243–0.260] | 0.151 | 0.38 / 0.33 |
| G2.2 | 0.455 [0.423–0.496] | 0.272 | 0.326 [0.304–0.353] | 0.192 | 0.42 / 0.36 |
| G3.1 | 0.382 [0.367–0.401] | 0.227 | 0.264 [0.256–0.276] | 0.161 | 0.42 / 0.35 |
| G3.2 | 0.448 [0.439–0.458] | 0.259 | 0.314 [0.305–0.321] | 0.184 | 0.44 / 0.35 |
| all 216 | 0.427 [0.342–0.496] | 0.231 | 0.305 [0.243–0.353] | 0.165 | — |

Per-halo ranges are 0.01–0.08 wide; halo medians span 0.35–0.46. **The six
halos are six correlated clusters, not 216 draws** — the same thing the
Stage 4b continuum test found (leave-one-halo-out fails; effective N = 3
Local Group realisations). G2.2's stars (0.272; particle file G2.1) are slightly *more*
concentrated than F98/Cao (0.267/0.269): "close to the bulge templates" is
the right phrase for the stellar side, not "less concentrated than".

## Suggested SM table (copy-paste; compiles in revtex4-2, single column)

Test-compiled: the first draft was 17 pt too wide for one column; this
version uses primed headers (defined in the caption), short template
names, and `\footnotesize`, and compiles with no overfull box. Without
`\footnotesize` it is still 6 pt over; `table*` (full width) also fits if
the smaller font is unwanted.

```latex
\begin{table}[b]
\caption{Concentration $C\equiv F(<3^{\circ})/F(<10^{\circ})$ at $2^{\circ}$
Gaussian smoothing for the annihilation ($J$) and old-stellar ($\star$) maps
of each halo: the median over its 36 observer azimuths, with the full range
for $J$ in parentheses. The right-hand columns excise the inner
$1.6^{\circ}$, the angular scale of the softening length, from both
numerator and denominator:
$C'\equiv F(1.6^{\circ}<\theta<3^{\circ})/F(1.6^{\circ}<\theta<10^{\circ})$.
The bulge templates and the spherical gNFW$^{2}$ profile are listed for
comparison.}
\label{tab:conc}
\footnotesize
\begin{ruledtabular}
\begin{tabular}{lcccc}
 & $C_{J}$ & $C_{\star}$ & $C'_{J}$ & $C'_{\star}$ \\
\hline
G1.1 & 0.46 (0.44--0.47) & 0.22 & 0.33 (0.32--0.34) & 0.16 \\
G1.2 & 0.38 (0.36--0.44) & 0.22 & 0.27 (0.25--0.31) & 0.16 \\
G2.1 & 0.35 (0.34--0.36) & 0.21 & 0.25 (0.24--0.26) & 0.15 \\
G2.2 & 0.46 (0.42--0.50) & 0.27 & 0.33 (0.30--0.35) & 0.19 \\
G3.1 & 0.38 (0.37--0.40) & 0.23 & 0.26 (0.26--0.28) & 0.16 \\
G3.2 & 0.45 (0.44--0.46) & 0.26 & 0.31 (0.30--0.32) & 0.18 \\
\hline
All 216 views & 0.43 (0.34--0.50) & 0.23 & 0.31 (0.24--0.35) & 0.17 \\
\hline
Freudenreich & 0.27 & & 0.20 & \\
Cao & 0.27 & & 0.19 & \\
Coleman & 0.29 & & 0.21 & \\
gNFW$^{2}$ ($\gamma=1.2$) & 0.51 & & 0.30 & \\
\end{tabular}
\end{ruledtabular}
\end{table}
```

## Suggested SM text (after the concentration paragraph in §A.2)

```latex
Table~\ref{tab:conc} lists the concentrations by halo, together with $C'$, the same ratio with the inner $1.6^{\circ}$ excised. The spread among the
36 azimuths of a single halo ($0.01$--$0.08$) is smaller than the spread
among halos (medians $0.35$--$0.46$), so the 216 views are six correlated
samples rather than 216 independent ones. Between $38\%$ and $44\%$ of
$F(<3^{\circ})$ in the annihilation maps arises within $1.6^{\circ}$ of the
center, the angular scale of the softening length, and $33$--$36\%$ for the
stellar maps; these are upper bounds on the unresolved contribution, since
the same lines of sight also traverse well-resolved material at larger
radii. Excising this region from both numerator and denominator lowers
every concentration but changes no ranking: every annihilation map remains
more concentrated than each bulge template and than its own stellar map
(216 of 216) for exclusion radii of $1^{\circ}$, $1.6^{\circ}$, and
$2^{\circ}$, at both $2^{\circ}$ and $3^{\circ}$ smoothing. The conclusion
therefore does not rest on the unresolved center.
```

Optional closing sentence (smoothing-dependent, see the warning above):

```latex
With the center excised, the simulated annihilation maps are as
concentrated as the spherical gNFW$^{2}$ profile at $2^{\circ}$ smoothing
(118 of 216 views exceed it) but not at $3^{\circ}$ (9 of 216).
```

The main-text sentence "smoothing does not remove this unresolved
contribution to the central flux" can now be followed by a pointer: "; the
Supplemental Material shows that excising it does not change the result."
