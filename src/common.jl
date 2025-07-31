module Common

"""
    canonical(x)

Return a canonical form of `x`.

For some types, there is more than value of the type that represents a single semantic value.
In this case, `canonical(x)` maps, in general, multiple type values of `x` with the same
semantic value to a single canonical type value.

Note that `Rational(0, a)` is canonicalized to `Rational(0, 1)` on construction.
"""
canonical(x) = x

function one_over_root_two(::Type{T}) where {T <: Real}
    one(T) / sqrt(T(2))
end

"""
    imaginary(::Type{T}) where {T <: Real}
    imaginary(::Type{Complex{T}}) where T

The imaginary unit as type `T`.

If the imaginary unit cannot be represented as `T`, for example
when `T` is `Real`, then `Complex{T}` is returned.
"""
function imaginary(::Type{T}) where {T <: Real}
    complex(zero(T), one(T))
end

function imaginary(::Type{Complex{T}}) where T
    imaginary(T)
end

"""
    sqrt_imaginary(::Type{T})

The principal square root of the imaginary unit as type `T`.

If the imaginary unit cannot be represented as `T`, for example
when `T` is `Real`, then `Complex{T}` is returned.
"""
function sqrt_imaginary(::Type{T}) where {T <: Real}
    cispi(one(T)/4)
end

function sqrt_imaginary(::Type{Complex{T}}) where {T <: Real}
    cispi(one(T)/4)
end

function sqrt_imaginary(::Type{BigFloat})
    cispi(one(BigFloat)/4)
end

end # module Common
