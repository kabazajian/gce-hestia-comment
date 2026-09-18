# Stage 3a anchor: run Muru et al.'s OWN projection code (gce_in_hestia,
# b50114d9, unmodified sources) on the G3.2 particle file, with the hard-coded
# Appendix A parameters, so our Python `project_muru` can be validated against
# their pipeline numerically. See validation/muru_redo_plan_2026-08-30.md §3a.
#
# Geometry is supplied explicitly (their AHF loaders are private; the stub
# HestiaUtils errors if touched). The frame vectors are deliberately
# NON-axis-aligned so a handedness or rotation error cannot hide.
#
# Run:  julia --project=. -t 4 anchor_muru_g32.jl

using Arrow, DataFrames, DelimitedFiles, LinearAlgebra

const THEIR_SRC = get(ENV, "GCE_IN_HESTIA", "/Users/aba/work/gce_in_hestia") * "/src"
include(joinpath(THEIR_SRC, "GalacticCoordinates.jl"))   # verbatim
include(joinpath(THEIR_SRC, "BulgeDensityProjection.jl")) # verbatim
using .GalacticCoordinates
using .BulgeDensityProjection

const H = 0.677
const OUT = joinpath(@__DIR__, "out")
mkpath(OUT)

# --- geometry (identical constants in anchor_compare.py — keep in sync) -----
# G3.2 AHF center, absolute box Mpc/h (Muru email 2026-08-30), -> kpc/h to
# match the kpc/h working units their coordsinMpc=true conversion produces.
p_c = [46.7533242, 50.3214297, 47.7923984] .* 1e3
# Deliberately non-axis-aligned, deterministic:
v_n = normalize([0.2, 0.3, 0.9])
v_d_raw = [1.0, 0.1, -0.05]
v_d = normalize(v_d_raw .- (v_d_raw ⋅ v_n) .* v_n)   # in-plane, perp to v_n
const ANGLE = 0.0                    # their hard-coded figure-path azimuth
const D_GC = 8.0 * H                 # kpc/h — their hard-coded 8 kpc
const D_MAX = 15.0 * H               # kpc/h — their hard-coded 15 kpc cut
const APERTURE = 3.0                 # deg — their hard-coded search distance
# Anchor grid: their ±20° figure grid at their exact 0.1°
const GRID = (-20.0, 20.0, 0.1)

t0 = time()
particles = Arrow.Table("/Users/aba/work/hestia_muru/HESTIA_37_11_8192_G32_minimal_particles.arrow") |> DataFrame |> copy
println("read $(nrow(particles)) rows  [$(round(time()-t0, digits=1)) s]")

# their dm filter (calc_methods_densityprojections.jl:19); ptype may
# deserialize as String or Symbol depending on Arrow round-trip — accept both,
# and require the count to be the independently-measured 4,472,132.
filter!(:ptype => p -> Symbol(p) == :dm, particles)
@assert nrow(particles) == 4_472_132 "dm count mismatch: $(nrow(particles))"

# their figure path: coords -> filter d (calc_density_maps lines 58-61)
add_galacticcoordinates!(particles, v_n, v_d, p_c, ANGLE, D_GC)
p_obs = generate_sun_position(p_c, v_n, v_d, ANGLE, D_GC)  # same α ⇒ same point
filter!(:d => <=(D_MAX), particles)
println("after 15 kpc cut: $(nrow(particles)) particles  [$(round(time()-t0, digits=1)) s]")

lbs_grid = generate_sym_lbsgrid(GRID...)
dens = calc_density_on_grid(particles, lbs_grid, APERTURE)
println("map done: $(size(dens))  [$(round(time()-t0, digits=1)) s]")

writedlm(joinpath(OUT, "julia_losmass.tsv"), dens)
writedlm(joinpath(OUT, "julia_observer.tsv"), p_obs')
# first 200k particles' (l,b,d) for a transform-level cross-check
n = min(200_000, nrow(particles))
writedlm(joinpath(OUT, "julia_lbd_head.tsv"),
         [particles.ParticleIDs[1:n] particles.l[1:n] particles.b[1:n] particles.d[1:n]])
println("wrote $(OUT)  [$(round(time()-t0, digits=1)) s]")
