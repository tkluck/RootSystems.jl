RootSystems.jl
==============

A simple library that contains implementations of the [irreducible root systems](https://en.wikipedia.org/wiki/Root_system).

Example use:

```julia
julia> using RootSystems

julia> Φ = D(4)
D{4}()

julia> length(dynkin_diagram_automorphisms(Φ))
6

julia> α, β = roots(Φ); # first two roots in iteration order

julia> coordinates(α)
4-element Array{Int64,1}:
 1
 1
 0
 0

julia> coordinates(β)
4-element Array{Int64,1}:
  1
 -1
  0
  0

julia> reflectalong(coordinates(α), coordinates(β))
4-element Array{Int64,1}:
  1
 -1
  0
  0

```

This library was created for a specific computation [related to this
paper](https://aip.scitation.org/doi/abs/10.1063/1.4705269) and is therefore
not as general as it could/should be.
