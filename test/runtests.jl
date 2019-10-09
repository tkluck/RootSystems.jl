using Test
using RootSystems

const TestSystems = [
    [A(n) for n in 1:10];
    [D(n) for n in 1:10];
    [E(n) for n in 6:8];
]

@testset "RootSystems.jl" begin
    @testset "Consistency in $Φ" for Φ in TestSystems
        @test length(Φ) == length(roots(Φ))
        @test length(Φ) == 2length(positive_roots(Φ))
        @test rank(Φ)   == length(simple_roots(Φ))

        @test Set(positive_roots(Φ)) == Set(filter(β -> all(coefficients_on_simple_roots(Φ, β) .>= 0), roots(Φ)))

        Δ = map(collect∘coordinates, simple_roots(Φ))
        @test all(roots(Φ)) do β
            sum(c*δ for (c,δ) in zip(coefficients_on_simple_roots(Φ, β), Δ)) == collect(coordinates(β))
        end
    end
end
