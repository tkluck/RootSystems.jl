using LinearAlgebra: I

struct E6 <: RootSystem end
struct E7 <: RootSystem end
struct E8 <: RootSystem end

const E = Union{E6, E7, E8}

E(n) = n == 6 ? E6() : n == 7 ? E7() : n == 8 ? E8() : error("There is no E$n root system")


struct HalfIntegerE8Root <: Root
    signs :: UInt8
end

E8Root = Union{DRoot{8}, HalfIntegerE8Root}

Base.:-(α::HalfIntegerE8Root) = typeof(α)(~α.signs)

rank(::E6) = 6
rank(::E7) = 7
rank(::E8) = 8

Base.length(::E6) =  72
Base.length(::E7) = 126
Base.length(::E8) = 240

coordinates(α::HalfIntegerE8Root) = [
    !iszero(α.signs & (256 >> k)) ? -1//2 : 1//2
    for k in 1:8
]

function simple_roots(Φ::E8)
    return E8Root[
        simple_roots(D(8))[1:6];
        DRoot{8}(6, 7, +1, +1);
        HalfIntegerE8Root(255);
    ]
end

simple_roots(Φ::E7) = simple_roots(E(8))[2:end]
simple_roots(Φ::E6) = simple_roots(E(8))[3:end]

function roots(Φ::E8)
    return E8Root[
        roots(D(8));
        [
            HalfIntegerE8Root(c)
            for c in 0:255
            if iseven(count_ones(c))
        ]
    ]
end

roots(Φ::E7) = filter(r -> iszero(coefficients_on_simple_roots(E(8), r)[1]), roots(E(8)))
roots(Φ::E6) = filter(r -> iszero(coefficients_on_simple_roots(E(7), r)[1]), roots(E(7)))

function coefficients_on_simple_roots(Φ::E8, α::DRoot{n}) where n
    D_coeffs = coefficients_on_simple_roots(D(8), α)

    return Int[D_coeffs[1:6]; 0; 0] .+
        D_coeffs[7] .* [ 1, 2, 3, 4, 5, 2, 4, 2] .+
        D_coeffs[8] .* [-1,-2,-3,-4,-5,-3,-3,-2]
end

function coefficients_on_simple_roots(Φ::E8, α::HalfIntegerE8Root) where n
    num_D_roots = count_zeros(α.signs) ÷ 2
    c = α.signs
    k = 0
    D_roots = map(1:num_D_roots) do _
        δ = leading_ones(c) + 1
        i = (k += δ)
        c <<= δ
        δ = leading_ones(c) + 1
        j = (k += δ)
        c <<= δ
        DRoot{8}(i, j, +1, +1)
    end

    return mapreduce(β -> coefficients_on_simple_roots(Φ, β), +, D_roots, init=[0,0,0,0,0,0,0,1]);
end

function coefficients_on_simple_roots(Φ::E7, α)
    coeff = coefficients_on_simple_roots(E(8), α)
    iszero(coeff[1]) || error("$α is not in $Φ")
    return coeff[2:end]
end

function coefficients_on_simple_roots(Φ::E6, α)
    coeff = coefficients_on_simple_roots(E(7), α)
    iszero(coeff[1]) || error("$α is not in $Φ")
    return coeff[2:end]
end

struct E6Identity end
struct E7Identity end
struct E8Identity end
struct E6Flip end

Base.:∘(::Id, ::Id) where Id <: Union{E6Identity, E7Identity, E8Identity} = Id()
Base.:∘(::E6Identity, ::E6Flip) = E6Flip()
Base.:∘(::E6Flip, ::E6Identity) = E6Flip()

dynkin_diagram_automorphisms(::E6) = [
    E6Identity(),
    E6Flip(),
]

dynkin_diagram_automorphisms(::E7) = [EIdentity()]
dynkin_diagram_automorphisms(::E8) = [EIdentity()]

coordinates(::E6Identity) = I
coordinates(::E7Identity) = I
coordinates(::E8Identity) = I
coordinates(::E6Flip) = [
     1  1 -3 0  0  0  0  1
     1  1 -3 0  0  0  0  1
    -3 -3 -3 0  0  0  0 -3
     0  0  0 3  3  3  3  0
     0  0  0 3  3 -3 -3  0
     0  0  0 3 -3  3 -3  0
     0  0  0 3 -3 -3  3  0
     1  1 -3 0  0  0  0  1
]//6
