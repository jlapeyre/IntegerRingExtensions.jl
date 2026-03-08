"""
module Angles2

Exports:
$(EXPORTS)
"""
module Angles2

using DocStringExtensions: EXPORTS, TYPEDFIELDS, SIGNATURES, TYPEDSIGNATURES, TYPEDEF
import Random

export Dar, radtodar, scalepi, unscalepi, minus_one_to_one, random_angle, random_phase

const PRETTY = MIME"text/plain"

"""
    $(TYPEDEF)

An angular unit equal to radians scaled by `π`.

# Fields
$(TYPEDFIELDS)

# Examples
```jldoctest
julia> Dar(3)
3 π

julia> Dar(3) == 3pi
true
```
`Dar(1)` is exactly equal to `π`
```jldoctest
julia> 3 * Dar(1//3)
1//1 π

julia> float(pi) == pi
false

julia> Dar(1) == pi
true

julia> 3 * Dar(1//3)
1//1 π

julia> 3 * Dar(1//3) == pi
true
```
"""
struct Dar{T} #  <: Real I think I don't want this
    """coefficient of π"""
    c::T

    function Dar(x)
        new{typeof(x)}(x)
    end
end

coeff(d::Dar) = d.c

function Base.show(io::IO, ::PRETTY, d::Dar)
    print(io, coeff(d), " π")
end

Base.AbstractFloat(::Type{Dar{T}}) where {T} = float(T)
Base.AbstractFloat(a::Dar{T}) where {T} = unscalepi(float(coeff(a)))
Base.Float64(a::Dar) = unscalepi(float(coeff(a)))

# Base.convert(::Type{Float64}, a::Dar) = Float64(a)

# These are odd. Don't seem to make sense.
# Return type is different from input type.
# But, sometimes, we have `angle::T` and want `one(T)`.
# This allows using both `angle::Float64` and `Dar{T}` in one code path.
Base.one(::Type{Dar{T}}) where {T} = one(T)
Base.one(::Dar{T}) where {T} = one(T)

# This is sensible
Base.zero(::Type{Dar{T}}) where {T} = Dar(zero(T))
Base.zero(::Dar{T}) where {T} = Dar(zero(T))

Base.:+(d::Dar, x::Real) = Dar(d.c + scalepi(x))
Base.:-(d::Dar, x::Real) = Dar(d.c - scalepi(x))
Base.:+(x::Real, d::Dar) = d + x
Base.:-(x::Real, d::Dar) = d - x
Base.:+(d1::Dar, d2::Dar) = Dar(d1.c + d2.c)
Base.:-(d1::Dar, d2::Dar) = Dar(d1.c - d2.c)
Base.:-(a::Dar) = Dar(-a.c)

Base.:(==)(a::Dar, b::Dar) = a.c == b.c

Base.:(==)(a::Dar, x::Real) = unscalepi(a.c) == x
Base.:(==)(x::Real, a::Dar) = x == a

# float(pi) != pi
# But we test if a is exactly pi.
Base.:(==)(a::Dar{T}, ::Irrational{:π}) where {T<:Union{Integer, Rational}} = isone(a.c)

# Note that we don't want ::Dar * ::Dar.
# That would make no sense
Base.:*(x::Real, d::Dar) = Dar(x * d.c)
Base.:*(d::Dar, x::Real) = x * d

Base.:/(d::Dar, x::Real) = Dar(d.c / x)
Base.://(d::Dar, x::Real) = Dar(d.c // x)
Base.div(d::Dar, x::Real) = Dar(div(d.c, x))
Base.:/(d::Dar, ::Irrational{:π}) = d.c

Dar(d::Dar) = d

## We might want to use these rotate the angle into a given range.
## For example, some math proofs depend on this.

"""
    $(SIGNATURES)

Add an even integer to `x` such that the result is in `[-1, 1]` and return the result.
"""
function minus_one_to_one(x)
    (f, w) = modf(x) # fractional, whole
    iszero(f) && return iseven(x) ? zero(x) : one(x)
    isodd(w) && (w += sign(x))
    return x - w
end

"""
    $(SIGNATURES)

Add an even integer to `x` such that the result is in `[0, 2]` and return the result.
"""
function zero_to_two(x)
    mod(x, 2)
end

# Probably want a Base.convert
"""
    radtodar(theta)::Dar

Convert angle `theta` (in radians) to `Dar` of same value.

# Examples
```jldoctest
julia> d = radtodar(3.0)
0.954929658551372 π

julia> float(d) == 3.0
true
```
"""
function radtodar(theta)
    Dar(theta / pi)
end

scalepi(x::Real) = x / pi
scalepi(x::Dar) = x
unscalepi(x::Real) = x * pi
unscalepi(x::Dar) = x

# intdiv(x::Integer, n::Integer) = div(x, n)
# intdiv(x::Real, n::Integer) = x / n
# #intdiv(x::Rational, n::Integer) = x * 1//n
# intdiv(x::Dar, n::Integer) = Dar(intdiv(x.c, n))

cscpi(x) = inv(sinpi(x))
secpi(x) = inv(cospi(x))

import Base: cos, sin, cis, tan, sincos, csc, sec

for func in (:cos, :sin, :cis, :tan, :sincos, :csc, :sec)
    funcpi = Symbol(func, :pi)
    @eval $func(dar::Dar) = $funcpi(dar.c)
end

# A random angle between -pi and pi
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar{T}}) where {T<:AbstractFloat}
    Dar(2*rand(T) - 1)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar})
    T = Float64
    Dar(2*rand(T) - 1)
end

Base.isapprox(a::Dar, b::Dar; kws...) = isapprox(a.c, b.c; kws...)
Base.isapprox(a::Dar, x::Real; kws...) = isapprox(unscalepi(a.c), x; kws...)
Base.isapprox(x::Real, a::Dar ; kws...) = isapprox(a, x; kws...)


"""
    random_angle(::Type{T}=Float64, shape...)::T where {T}

Sample from the uniform distribution on `[-pi, pi]`
"""
random_angle(args...) = random_angle(Float64, args...)

function random_angle(::Type{T}, args...) where {T<:AbstractFloat}
    flts = rand(T, args...)
    f = x -> 2 * T(pi) * (x - T(1)/2)
    if isempty(args)
        return f(flts)
    end
    map!(f, flts, flts)
end

function random_angle(::Type{T}, args...) where {T<:Dar}
    rand(T, args...)
end

"""
    random_phase(::Type{T}=Float64) where {T}

Return a random phase, that is `exp(i π α)` where α is a random angle.
"""
function random_phase(::Type{T}=Float64) where {T}
    cis(random_angle(T))
end

end # module Angles2
