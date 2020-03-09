struct D{n} <: RootSystem end

D(n) = D{n}()

struct DRoot{n} <: Root
    i :: Int
    j :: Int
    sign_i :: Int
    sign_j :: Int
end

Base.:-(α::DRoot) = typeof(α)(α.i, α.j, -α.sign_i, -α.sign_j)

rank(::D{n}) where n = n
Base.length(::D{n}) where n = 2n*(n-1)

coordinates(α::DRoot{n}) where n = [
    k == α.i ? α.sign_i : k == α.j ? α.sign_j : 0
    for k in 1:n
]

function simple_roots(Φ::D{n}) where n
    return [map(i->DRoot{n}(i, i+1, +1, -1), 1:(n-1)); DRoot{n}(n-1, n, +1, +1)]
end

function positive_roots(Φ::D{n}) where n
    return [
        [
            DRoot{n}(i, j, +1, -1)
            for i in 1:n
            for j in i+1:n
        ]; [
            DRoot{n}(i, j, +1, +1)
            for i in 1:n
            for j in i+1:n
        ]
    ]
end

function roots(Φ::D{n}) where n
    return [
        DRoot{n}(i, j, sign_i, sign_j)
        for i in 1:n
        for j in i+1:n
        for sign_i in (+1, -1)
        for sign_j in (+1, -1)
    ]
end

function coefficients_on_simple_roots(Φ::D{n}, α::DRoot{n}) where n
    if (α.sign_i, α.sign_j) == (+1, -1)
        res = zeros(Int, n)
        for k in 1:n-1
            if α.i <= k < α.j
                res[k] = 1
            end
        end
        return res
    elseif (α.sign_i, α.sign_j) == (+1, +1)
        res = zeros(Int, n)
        if α.j == n
            res[α.i:n-2] .= 1
            res[n] = 1
        else
            for k in 1:n
                if α.i <= k < α.j
                    res[k] = 1
                elseif α.j <= k < n-1
                    res[k] = 2
                elseif n-1 <= k <= n
                    res[k] = 1
                end
            end
        end
        return res
    else
        return -coefficients_on_simple_roots(Φ, -α)
    end
end

struct DAut{n}
    permutation :: Vector{Int}
    signs       :: Vector{Int}
end

function Base.:∘(f::Aut, g::Aut) where Aut <: DAut
    Aut(
        f.permutation[g.permutation],
        f.signs[g.permutation] .* g.signs
    )
end

Base.:(==)(a::Aut, b::Aut) where Aut <: DAut = a.permutation == b.permutation && a.signs == b.signs

weyl_group(::D{n}) where n = (
    DAut{n}(p, [iszero(s & big"1" << b) ? +1 : -1 for b in 0 : n-1])
    for p in permutations(1:n)
    for s in 0:(big"2"^n - 1)
    if iseven(count_ones(s))
)

dynkin_diagram_automorphisms(::D{n}) where n = [
    DAut{n}([1:n;], [[+1 for _ in 1:n-1]; +1]),
    DAut{n}([1:n;], [[+1 for _ in 1:n-1]; -1]),
]

function coordinates(f::DAut{n}) where n
    res = spzeros(Int, n, n)
    for (i, (j, s)) in enumerate(zip(f.permutation, f.signs))
        res[j, i] = s
    end
    return res
end

struct D4Aut
    permutation :: Vector{Int}
end

Base.:(==)(a::D4Aut, b::D4Aut) = a.permutation == b.permutation

Base.:∘(f::D4Aut, g::D4Aut) = D4Aut(f.permutation[g.permutation])

dynkin_diagram_automorphisms(::D{4}) = [
    D4Aut([1, 2, 3, 4]),
    D4Aut([1, 2, 4, 3]),
    D4Aut([3, 2, 1, 4]),
    D4Aut([3, 2, 4, 1]),
    D4Aut([4, 2, 1, 3]),
    D4Aut([4, 2, 3, 1]),
]

function coordinates(f::D4Aut)
    Δ = map(coordinates, simple_roots(D(4)))
    id = hcat(Δ...)' // 1
    x = hcat(Δ[f.permutation]...)' // 1
    return x \ id
end
