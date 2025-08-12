module Common

# TODO: cleanup, prolly remove, root_two, etc. in favor of T(singleton)
# Design using `root_two` agrees Julia Base use of zero() and one().
# But, I think using T(singleton) is more robust and flexible.
# It allows defining method for, say :*(::SingletonType, obj)
# The conversions are decide on in the :* method, where they belong.

export canonical, isrational, coeffs, params, isunit, invstrict

"""
    canonical(x)

Return a canonical form of `x`.

For some types, there is more than value of the type that represents a single semantic value.
In this case, `canonical(x)` maps, in general, multiple type values of `x` with the same
semantic value to a single canonical type value.

Note that `Rational(0, a)` is canonicalized to `Rational(0, 1)` on construction.
"""
canonical(x) = x

canonical(a::AbstractArray) = map(canonical, a)

"""
    one_over_root_two(::Type{T}) where {T <: Number}

Return a value of type `T` representing the reciprocal of the square root of two.
"""
function one_over_root_two(::Type{T}) where {T <: Number}
    one(T) / sqrt(T(2))
end

"""
    root_two(::Type{T}) where {T <: Number}

Return a value of type `T` representing the square root of two.
"""
function root_two(::Type{T}) where {T <: Number}
    sqrt(T(2))
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

# Get coefficients of instance of type
function coeffs end

# Get parameters of instance of type
function params end

##
## Following are fallback methods for a few arithmetic operations.
## Specialized methods exist for ring types.
##

"""
    mul_root_two(x)

Return `x` multiplied by the square root of two.
"""
function mul_root_two(x)
    x * sqrt(typeof(x)(2))
end

"""
    mul_one_over_root_two(x)

Return `x` multiplied by the reciprocal of the square root of two.

For some rings, all elements may be divided by the square root of two.
However, we call this multiplication to emphasize that this function can be called on elements of a ring.
"""
function mul_one_over_root_two(x)
    x * inv(sqrt(typeof(x)(2)))
end


"""
    mul_half(x, n::Integer=1)

Return `x` multiplied by the reciprocal of two to the `n`th power.
"""
function mul_half(x, n::Integer=1)
    x * (inv(typeof(x)(2))^n)
end

function mul_two(x, n::Integer=1)
    x * (typeof(x)(2))^n
end

"""
    conj_root_two(x::Number)

Represents the automorphism mapping `√2` to `-√2` for rings containing `√2`.

For elements of other rings, `conj_root_two` is the identity.
"""
conj_root_two(x::Number) = x

function norm_root_two end

"""
    isrational(x)

Return `true` iff the value `x` represents a rational number or a gaussian rational number.
"""
function isrational end

isrational(x::Number) = true
isrational(x::AbstractIrrational) = false

"""
   isunit(x::T)::Bool where {T}

Return `true` if there exists a `y::T` such that `x * y == one(T)`.
"""
function isunit end

isunit(r::Rational) = !iszero(r)
isunit(x::AbstractFloat) = !iszero(x)
isunit(c::Complex{T}) where {T<:AbstractFloat} = !iszero(c)
isunit(n::Integer) = isone(n) || isone(-n)

"""
    invstrict(x::T)

Return `y::T` such that `x * y == one(T)`.

If no such `y` exists, throw an error.
"""
function invstrict end

invstrict(x::Union{AbstractFloat}) = inv(x)
invstrict(x::Complex{T}) where {T <: Union{AbstractFloat, Rational}} = inv(x)

# Base.inv(::Rational) returns 1//0 for inv(0//n)
function invstrict(r::Rational)
    iszero(r) && throw(ArgumentError(lazy"$r has no inverse of type $(typeof(r))"))
    inv(r)
end

function invstrict(n::Integer)
    isone(n) && return n
    isone(-n) && return n
    throw(ArgumentError(lazy"$n has no inverse of type $(typeof(n))"))
end

function invstrict(c::Complex{T}) where {T<:Integer}
    isone(c) && return c
    isone(-c) && return c
    (r, i) = reim(c)
    if iszero(r)
        isone(i) && return Complex{T}(0, -one(T))
        isone(-i) && return Complex{T}(0, one(T))
    end
    throw(ArgumentError(lazy"$c has no inverse of type $(typeof(c))"))
end

end # module Common
