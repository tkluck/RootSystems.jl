using Combinatorics: permutations
using SparseArrays: spzeros

struct A{n} <: RootSystem end

A(n) = A{n}()

struct ARoot{n} <: Root
    i :: Int
    j :: Int
end

rank(::A{n}) where n = n
Base.length(::A{n}) where n = (n+1)*n

coordinates(α::ARoot{n}) where n = [
    k == α.i ? 1 : k == α.j ? -1 : 0
    for k in 1:n+1
]

function simple_roots(Φ::A{n}) where n
    return map(i->ARoot{n}(i,i+1), 1:n)
end

function positive_roots(Φ::A{n}) where n
    return [
        ARoot{n}(i, j)
        for i in 1:n
        for j in i+1:n+1
    ]
end

function roots(Φ::A{n}) where n
    return [
        ARoot{n}(i, j)
        for i in 1:n+1
        for j in 1:n+1
        if i != j
    ]
end

function coefficients_on_simple_roots(Φ::A{n}, α::ARoot{n}) where n
    res = zeros(Int, n)
    for k in 1:n
        i_seen = α.i <= k
        j_seen = α.j <= k

        if i_seen && !j_seen
            res[k] = 1
        elseif !i_seen && j_seen
            res[k] = -1
        end
    end
    return res
end

struct AAut{n}
    permutation :: Vector{Int}
    sign        :: Int
end

function Base.:∘(f::Aut, g::Aut) where Aut <: AAut
    Aut(
        f.permutation[g.permutation],
        f.sign .* g.sign,
    )
end

weyl_group(::A{n}) where n = (
    AAut{n}(p, +1)
    for p in permutations(1:n + 1)
)

dynkin_diagram_automorphisms(::A{n}) where n = [
    AAut{n}([1:n + 1;],    +1),
    AAut{n}([n + 1:-1:1;], -1),
]

function coordinates(f::AAut{n}) where n
    res = spzeros(Int, n+1, n+1)
    for (i, j) in enumerate(f.permutation)
        res[j, i] = f.sign
    end
    return res
end
