```@meta
CurrentModule = RootOnes
```

# RootOnes

Package [RootOnes](https://github.com/jlapeyre/IntegerRingExtensions.jl/tree/main/lib/RootOnes) providing roots of unity.

```@docs
RootOnes
```

### Types
```@docs
RootOne
Omega
omega
```

### Methods
```@docs
isprimitive
```

### Methods for imported functions
```@docs
imaginary(::Type{RootOne{D}}) where {D}
sqrt_imaginary(::Type{RootOne{D}}) where {D}
isunit(::RootOne)
conj_root_two(r::RootOne{8})
norm_root_two(r::RootOne{8})
Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{RootOne{N}}) where {N}
Base.inv(::RootOne)
```
