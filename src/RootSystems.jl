module RootSystems

abstract type RootSystem end
abstract type Root end
abstract type Automorphism end

function rank end
function coordinates end
function simple_roots end
function positive_roots end
function roots end
function coefficients_on_simple_roots end

function weyl_group end
function dynkin_diagram_automorphisms end

function positive_roots(Φ::RootSystem)
    return filter(roots(Φ)) do β
        all(c -> c >= 0, coefficients_on_simple_roots(Φ, β))
    end
end

coordinates(a::AbstractVector) = a
⋅(a, b) = transpose(coordinates(a)) * coordinates(b)
reflectalong(α, β) = β - α⋅β*α

include("DynkinDiagramAutomorphisms.jl")
include("A-series.jl")
include("D-series.jl")
include("E-series.jl")

export RootSystem, A, D, E
export roots, simple_roots, positive_roots, coefficients_on_simple_roots, coordinates, rank
export dynkin_diagram_automorphisms
export ⋅, reflectalong


end # module
