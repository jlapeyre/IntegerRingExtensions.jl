# RootOnes.

[![Build Status](https://github.com/jlapeyre/RootOnes..jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jlapeyre/RootOnes..jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jlapeyre/RootOnes..jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jlapeyre/RootOnes..jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Provides `struct RootOne{N}`.


    struct RootOne{N}
    RootOne{N}(k::Integer)

`N`th roots of unity. `RootOne{N}(k)` is the `k`th power of the principal root.

`N` should always be a literal or a `const`. Otherwise performance is severely degraded.

`k` will be stored as `mod(k, N)`, which takes values in `(0,...,N-1)`.

The following alias is defined:
```julia
const Omega = RootOne{8}
```

### Examples

See also `Omega`.

```julia

julia> RootOne{8}(1)
RootOne{8}(1)
```

```julia
julia> RootOne{8}(-1)
RootOne{8}(7)

julia> -RootOne{8}(1) # Unary minus
RootOne{8}(5)

julia> -RootOne{8}(5)
RootOne{8}(1)

julia> -RootOne{9}(1)
ERROR: ArgumentError: InexactError unary minus of type RootOne{9}

julia> RootOne{8}(1)^8
RootOne{8}(0)

julia> isone(RootOne{8}(1)^8)
true

julia> one(RootOne{8})
RootOne{8}(0)

julia> one(RootOne{8}(5))
RootOne{8}(0)

julia> RootOne{8}(2) * RootOne{8}(3)
RootOne{8}(5)
```
