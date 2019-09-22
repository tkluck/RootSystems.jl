module RootSystems

abstract type RootSystem end
abstract type Root{System} end

function rank end
function coordinates end
function dynkin_diagram_automorphism_count end
function simple_roots end
function positive_roots end
function roots end
function coefficients_on_simple_roots end

function positive_roots(Φ::RootSystem)
    return filter(roots(Φ)) do β
        all(c -> c >= 0, coefficients_on_simple_roots(Φ, β))
    end
end

include("A-series.jl")
include("D-series.jl")
include("E-series.jl")

export RootSystem, A, D, E
export roots, simple_roots, positive_roots, coefficients_on_simple_roots, coordinates, rank
export dynkin_diagram_automorphism_count


end # module
