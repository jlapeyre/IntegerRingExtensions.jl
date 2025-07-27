```@meta
CurrentModule = RingExtensionsCommon
```

# RingExtensionsCommon

Documentation for [RingExtensionsCommon](https://github.com/jlapeyre/IntegerRingExtensions.jl/tree/main/lib/RingExtensionsCommon).

```@docs
RingExtensionsCommon
```

```@docs
canonical(x)
root_two(::Type{T}) where {T <: Number}
one_over_root_two(::Type{T}) where {T <: Number}
imaginary(::Type{T}) where {T <: Real}
sqrt_imaginary(::Type{T}) where {T <: Real}
mul_root_two(x)
mul_one_over_root_two(x)
mul_half(x, n::Integer=1)
mul_two(x, n::Integer=1)
conj_root_two(x::Number)
conj_root_D(x::Number)
norm_root_two
norm_root_D
invstrict
isrational
isunit
isimag
```
