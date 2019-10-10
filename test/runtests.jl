using Test
using RootSystems

using LinearAlgebra: I

const TestSystems = [
    [A(n) for n in 1:10];
    [D(n) for n in 2:10];
    [E(n) for n in 6:8];
]

@testset "RootSystems.jl" begin
    @testset "Consistency in $Φ" for Φ in TestSystems
        @test length(Φ) == length(roots(Φ))
        @test length(Φ) == 2length(positive_roots(Φ))
        @test rank(Φ)   == length(simple_roots(Φ))

        @test Set(positive_roots(Φ)) == Set(filter(β -> all(coefficients_on_simple_roots(Φ, β) .>= 0), roots(Φ)))

        Δ = map(coordinates, simple_roots(Φ))
        @test all(roots(Φ)) do β
            sum(c*δ for (c,δ) in zip(coefficients_on_simple_roots(Φ, β), Δ)) == coordinates(β)
        end

        if !(Φ isa E) && rank(Φ) < 4 # E not implemented; and restrict to small Weyl group
            W = collect(weyl_group(Φ))
            @test all(
                coordinates(f∘g) == coordinates(f)*coordinates(g)
                for f in rand(W, 100), g in rand(W, 100)
            )
        end

        if Φ != D(4) # not implemented
            G = dynkin_diagram_automorphisms(Φ)
            @test all( # g squares to identity
                coordinates(g) * coordinates(g) == I == coordinates(g ∘ g)
                for g in G
            )
            @test all( # g permutes the simple roots
                Set(coordinates(g) * δ for δ in Δ) == Set(Δ)
                for g in G
            )
        end
    end
end
