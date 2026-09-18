"""
STUB of Muru's private HestiaUtils (path = /z/moortis/.julia/dev/HestiaUtils in
their Manifest; not publicly available — CLAUDE.md rule 3). It exists ONLY so
that the unmodified `gce_in_hestia` sources load. Every function ERRORS if
called: the anchor drives their public compute path (`add_galacticcoordinates!`,
`generate_sun_position`, `calc_density_on_grid`) with explicitly supplied
centers and axes, and must never silently fall through to a fabricated loader.
Same UUID as their Manifest so `using HestiaUtils` resolves.
"""
module HestiaUtils

export SimulationSpecs, get_galcenter, get_galvectors_from_profile, read_ahfmergertree

struct SimulationSpecs
    simID
    snap
end

_stub(f) = error("HestiaUtils stub: `$f` is Muru's private data loader and is " *
                 "not available. The anchor must supply centers/axes explicitly.")

get_galcenter(args...) = _stub("get_galcenter")
get_galvectors_from_profile(args...) = _stub("get_galvectors_from_profile")
read_ahfmergertree(args...) = _stub("read_ahfmergertree")

end # module
