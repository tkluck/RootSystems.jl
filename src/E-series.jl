struct E6 <: RootSystem end
struct E7 <: RootSystem end
struct E8 <: RootSystem end

E(n) = n == 6 ? E6() : n == 7 ? E7() : n == 8 ? E8() : error("There is no E$n root system")


struct HalfIntegerE8Root
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

coordinates(α::HalfIntegerE8Root) = ntuple(8) do k
    !iszero(α.signs & (256 >> k)) ? -1//2 : 1//2
end

dynkin_diagram_automorphism_count(::E6) = 2
dynkin_diagram_automorphism_count(::E7) = 1
dynkin_diagram_automorphism_count(::E8) = 1

function simple_roots(Φ::E8)
    return E8Root[
        simple_roots(D(8))[1:6];
        DRoot{8}(6, 7, +1, +1);
        HalfIntegerE8Root(255);
    ]
end

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

