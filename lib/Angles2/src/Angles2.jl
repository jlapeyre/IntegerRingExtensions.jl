"""
    module Angles2

This module provides [`struct Dar`](@ref "Dar"), representing angles scaled by `π`. For example, `Dar(0) = Dar(2) == 1`

A couple of points:

* Rational multiples of `π` are represented exactly.
* `cos(Dar(x))` calls `cospi(x)`. Other trig functions are handled in the same way.

`cospi` can be more accurate. For example,
```julia-repl
julia> cosdiff(x) = cos(x * pi) - cospi(x);

julia> cosdiff(1/3)
1.1102230246251565e-16

julia> cosdiff(1.0000000001e10)
4.606315329169774e-13
```
So, if a function accepts `x::Union{Float64, Dar}` (or `x::Any`), then the caller has the option to use `cospi` wherever `cos` appears.

This is an experiment. How useful it is remains to be seen.

Method
- [`radtodar`](@ref)

Other functions
- [`minus_one_to_one`](@ref) -- shift input by integer such that it lies in `[-1, 1]`
- [`random_angle`](@ref) -- Generate random angles
- [`random_phase`](@ref) -- Generate random phase
"""
module Angles2

using DocStringExtensions: TYPEDFIELDS, SIGNATURES, TYPEDSIGNATURES, TYPEDEF
import Random

export Dar, radtodar, minus_one_to_one, random_angle, random_phase

const PRETTY = MIME"text/plain"

"""
    $(TYPEDEF)

An angular unit equal to radians scaled by `π`.

# Fields
$(TYPEDFIELDS)

### Examples
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

julia> Dar(3) / pi === 3
true

julia> div(Dar(6), 2)
3 π
```

A factor of π in the numerator and denominator cancels exactly.
```jldoctest
julia> Dar(4) / Dar(2)
2.0

julia> Dar(4) ÷ Dar(2)
2

julia> Dar(1) // Dar(2)
1//2
```

For improved accuracy, trigonometric functions call `cospi`, `sinpi`, etc. on the coeffient.
```jldoctest
julia> cos(Dar(1/2))
0.0

julia> cos(pi/2)
6.123233995736766e-17
```
"""
struct Dar{T<:Real} <: Number #  <: Real I think I don't want this
    """coefficient of π"""
    c::T

    function Dar(x::Real)
        new{typeof(x)}(x)
    end

    function Dar{T}(x::T) where T
        new{T}(x)
    end
end

Dar{T}(::Irrational{:π}) where {T<:Real} = Dar(one(T))
Dar(::Irrational{:π}) = Dar{Int}(1)
Dar{T}(d::Dar) where T = Dar{T}(T(d.c))
Dar(d::Dar) = d

coeff(d::Dar) = d.c

function Base.show(io::IO, ::PRETTY, d::Dar)
    print(io, coeff(d), " π")
end
Base.show(io::IO, d::Dar) = print(io, coeff(d), " π")

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

##
## TODO: Tame this promote_rule thicket
##

function Base.promote_rule(::Type{Dar{T}}, ::Type{S}) where {T<:Real,S<:Real}
    return Dar{promote_type(T, S)}
end

function Base.promote_rule(::Type{Dar{T}}, ::Type{Dar{S}}) where {T<:Real,S<:Real}
    return Dar{promote_type(T, S)}
end

function Base.promote_rule(::Type{Dar{T}}, ::Type{S}) where {T<:Real, S<:Irrational{:π}}
    return Dar{T}
end

function Base.promote_rule(::Type{S}, ::Type{Dar{T}}) where {T<:Real, S<:Irrational{:π}}
    return Dar{T}
end

function Base.promote_rule(::Type{Dar{T}}, ::Type{Irrational{:π}}) where {T<:Real}
    return Dar{T}
end

function Base.convert(::Type{T}, a::Dar) where {T<:AbstractFloat}
    return T(coeff(a)) * T(pi)
end

function Base.convert(::Type{Dar{T}}, x::AbstractFloat) where {T<:AbstractFloat}
    return Dar(T(x) / T(pi))
end

function Base.convert(::Type{Dar{T}},  a::V) where {T<:Real, V<:Irrational{:π}}
    return Dar(T(1))
end

Base.:-(a::Dar) = Dar(-a.c)
Base.:+(d1::Dar, d2::Dar) = Dar(+(promote(d1.c, d2.c)...))
Base.:-(d1::Dar, d2::Dar) = Dar(-(promote(d1.c, d2.c)...))
Base.:(==)(a::Dar, b::Dar) = a.c == b.c

Base.hash(d::Dar, h::UInt) = hash(d.c, h)
Base.isless(a::Dar, b::Dar) = isless(a.c, b.c)

# Note that we don't want ::Dar * ::Dar.
# That would make no sense
Base.:*(x::Real, d::Dar) = Dar(x * d.c)
Base.:*(d::Dar, x::Real) = x * d

Base.:/(d::Dar, x::Real) = Dar(d.c / x)
Base.://(d::Dar, x::Integer) = Dar(d.c // x)
Base.div(d::Dar, x::Integer) = Dar(div(d.c, x))
Base.:/(d::Dar, ::Irrational{:π}) = d.c

Base.:/(d1::Dar, d2::Dar) = d1.c / d2.c
Base.://(d1::Dar, d2::Dar) = d1.c // d2.c
Base.div(d1::Dar, d2::Dar) = div(d1.c, d2.c)

## We might want to use these rotate the angle into a given range.
## For example, some math proofs depend on this.

# We we are looking for agreement with
# myang(n) = angle(cispi(n))/pi
"""
    $(SIGNATURES)

Add an even integer to `x` such that the result is in `[-1, 1]` and return the result.

The choice at the endpoints is chosen to agree with `angle(cispi(x)) / pi`.

### Examples
```jldoctest
julia> minus_one_to_one(1.3)
-0.7000000000000002

julia> minus_one_to_one(2.3)
0.2999999999999998

julia> minus_one_to_one(-1.0)
-1.0

julia> minus_one_to_one(1.0)
1.0

julia> minus_one_to_one(-3.0)
-1.0
```
"""
function minus_one_to_one(x)
    if isinteger(x)
        if isodd(x)
            r = one(x)
        else
            r = zero(x)
        end
        if x < 0
            return -r
        else
            return r
        end
    end
    y = mod(x + 1, 2)
    return y - 1
end

"""
    $(SIGNATURES)

Add the unique even integer to `x` such that the result is in `[0, 2]` and return the result.
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

cotpi(x) = cospi(x) / sinpi(x)
Base.cot(d::Dar) = cotpi(d.c)

# A random angle between -pi and pi
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar{T}}) where {T<:AbstractFloat}
    Dar(T(2) * rand(rng, T) - T(1))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar})
    T = Float64
    Dar(T(2) * rand(rng, T) - T(1))
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
