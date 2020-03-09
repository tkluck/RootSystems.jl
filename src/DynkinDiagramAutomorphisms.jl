struct DynkinDiagramAutomorphism{RS}
    permutation :: Vector{Int}
end

function (f::DynkinDiagramAutomorphism{RS})(a) where RS <: RootSystem
    coeffs = coefficients_on_simple_roots(RS(), a)[f.permutation]
    return sum(c*r for (c,r) in zip(coeffs, simple_roots(RS())))
end
