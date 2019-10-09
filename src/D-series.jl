struct D{n} <: RootSystem end

D(n) = D{n}()

struct DRoot{n} #<: Root{D{n}}
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

dynkin_diagram_automorphism_count(::D{n}) where n = n == 4 ? 6 : 2

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

