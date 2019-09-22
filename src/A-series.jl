
struct A{n} <: RootSystem end

A(n) = A{n}()

struct ARoot{n} #<: Root{A{n}}
    i :: Int
    j :: Int
end

rank(::A{n}) where n = n
Base.length(::A{n}) where n = (n+1)*n

coordinates(α::ARoot{n}) where n = ntuple(Val(n+1)) do k
    k == α.i ? 1 : k == α.j ? -1 : 0
end

dynkin_diagram_automorphism_count(::A{n}) where n = 2

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

