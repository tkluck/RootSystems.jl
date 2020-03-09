using Test
using RootSystems

using LinearAlgebra: I

const TestSystems = [
    [A(n) for n in 1:10];
    [D(n) for n in 2:10];
    [E(n) for n in 6:8];
]

function op_power_by_squaring(op, x, n)
    n >= 1 || error("can't take operator power with n < 1")

    k = oftype(n, 2) ^ trailing_zeros(n)
    cur = x
    for _ in 2 : k
        cur = op(cur, cur)
    end
    res = cur

    while !iszero(k) && k <= n
        k <<= 1
        cur = op(cur, cur)

        if !iszero(n & k)
            res = op(res, cur)
        end
    end

    return res
end

function is_finite_group(op, G)
    length(G) == 0 && return true
    length(G) == 1 && op(first(G), first(G)) == first(G) && return true

    GG = Set(G)
    # closed
    all(op(a, b) in GG for a in G, b in G)
    # associativity
    all(op(a, op(b, c)) == op(op(a, b), c) for a in G, b in G, c in G) || return false
    # unit
    e = op_power_by_squaring(op, first(G), length(G))
    all(op(e, a) == op(a, e) == a for a in G) || return false
    # inverse
    all(op_power_by_squaring(op, a, length(G)) == e for a in G) || return false

    return true
end

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
            @test is_finite_group(∘, W)
        end

        G = dynkin_diagram_automorphisms(Φ)
        @test all( # g permutes the simple roots
            Set(coordinates(g) * δ for δ in Δ) == Set(Δ)
            for g in G
        )
        @test is_finite_group(∘, dynkin_diagram_automorphisms(Φ))
        if Φ != D(4)
            @test all( # g squares to identity
                coordinates(g) * coordinates(g) == I == coordinates(g ∘ g)
                for g in G
            )
        end
    end
end
