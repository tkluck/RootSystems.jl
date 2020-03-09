RootSystems.jl
==============

[![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url]


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

[travis-img]: https://travis-ci.org/tkluck/RootSystems.jl.svg?branch=master
[travis-url]: https://travis-ci.org/tkluck/RootSystems.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/mm4hem7506ct64ot?svg=true
[appveyor-url]: https://ci.appveyor.com/project/tkluck/rootsystems-jl
