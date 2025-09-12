module Angles

import Random
using ..Utils: PRETTY

"""
    Dar{T}

An angular unit that takes values on [-1, 1].
"""
struct Dar{T} #  <: Real I think I don't want this
    function Dar(y)
        x = minus_one_to_one(y)
        new{typeof(x)}(x)
    end
    x::T
end

function Base.show(io::IO, ::PRETTY, d::Dar)
    print(io, d.x, " π")
end

Base.AbstractFloat(::Type{Dar{T}}) where {T} = float(T)
Base.AbstractFloat(a::Dar{T}) where {T} = unscalepi(float(a.x))
Base.Float64(a::Dar) = unscalepi(float(a.x))

Base.convert(::Type{Float64}, a::Dar) = Float64(a)

Base.one(::Type{Dar{T}}) where {T} = one(T)
Base.one(::Dar{T}) where {T} = one(T)

Base.zero(::Type{Dar{T}}) where {T} = Dar(zero(T))
Base.zero(::Dar{T}) where {T} = Dar(zero(T))

Base.:+(d::Dar, x::Real) = Dar(d.x + scalepi(x))
Base.:-(d::Dar, x::Real) = Dar(d.x - scalepi(x))
Base.:+(x::Real, d::Dar) = d + x
Base.:-(x::Real, d::Dar) = d - x
Base.:+(d1::Dar, d2::Dar) = Dar(d1.x + d2.x)
Base.:-(d1::Dar, d2::Dar) = Dar(d1.x - d2.x)
Base.:-(a::Dar) = Dar(-a.x)
Base.:(==)(a::Dar, b::Dar) = a.x == b.x

Base.:(==)(a::Dar, x::Real) = unscalepi(a.x) == x
Base.:(==)(x::Real, a::Dar) = x == a

# Note that we don't want ::Dar * ::Dar.
# That would make no sense
Base.:*(x::Real, d::Dar) = Dar(x * d.x)
Base.:*(d::Dar, x::Real) = x * d

Dar(d::Dar) = d

"""
    minus_one_to_one(x)

Add an even integer to `x` such that the result is in `[-1, 1]` and return the result
"""
function minus_one_to_one(x)
    (f, w) = modf(x) # fractional, whole
    iszero(f) && return iseven(x) ? zero(x) : one(x)
    isodd(w) && (w += sign(x))
    return x - w
end

# Probably want a Base.convert
function radtodar(theta)
    Dar(theta / pi)
end

function radians(dar::Dar)
    dar.x * pi
end

scalepi(x::Real) = x / pi
scalepi(x::Dar) = x
unscalepi(x::Real) = x * pi
unscalepi(x::Dar) = x

intdiv(x::Real, n::Integer) = x / n
intdiv(x::Rational, n::Integer) = x * 1//n
intdiv(x::Dar, n::Integer) = Dar(intdiv(x.x, n))

cscpi(x) = inv(sinpi(x))
secpi(x) = inv(cospi(x))

import Base: cos, sin, cis, tan, sincos, csc, sec

for func in (:cos, :sin, :cis, :tan, :sincos, :csc, :sec)
    funcpi = Symbol(func, :pi)
    @eval $func(dar::Dar) = $funcpi(dar.x)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar{T}}) where {T<:AbstractFloat}
    Dar(2*rand(T) - 1)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Dar})
    T = Float64
    Dar(2*rand(T) - 1)
end

Base.isapprox(a::Dar, b::Dar; kws...) = isapprox(a.x, b.x; kws...)
Base.isapprox(a::Dar, x::Real; kws...) = isapprox(unscalepi(a.x), x; kws...)
Base.isapprox(x::Real, a::Dar ; kws...) = isapprox(a, x; kws...)


# Ang is nice. But Ang does not support BigFloat angles
# Because BigInt is not fixed precision.
# Dar above is more flexible

"""
    Ang{T<:Integer}

An angle in `(-π, π]`.

The angle is stored in field `x::T`.

# Examples
```jldoctest
julia> Ang(typemax(Int))
1.0 π

julia> Ang(typemin(Int))
-1.0 π
```
"""
struct Ang{T<:Integer}
    x::T
end

Ang(z::Float64) = Ang{Int}(z)

"""
    Ang{T}(z::Float64)::Ang{T} where {T<:Integer}

Return an `Ang` representing the angle `z π`.

`z` is shifted so that the result is in `(-π, π]`.

# Examples
```jldoctest
julia> Ang(0.5)
0.5 π

julia> cos(Ang(0.5))
0.0

julia> cos(101 * Ang(0.5))
0.0
```
"""
function Ang{T}(z::Float64) where {T<:Integer}
    zp = minus_one_to_one(z)
    if zp == 1
        Ang(typemax(T))
    elseif zp == -1
        Ang(-typemax(T))
    else
        Ang(round(T, zp * typemax(T)))
    end
end

function Base.:(==)(a::Ang{T},  b::Ang{T}) where {T}
    a.x == b.x
end

Base.:+(a::Ang{T}, b::Ang{T}) where {T} = Ang(a.x + b.x)
Base.:-(a::Ang{T}, b::Ang{T}) where {T} = Ang(a.x - b.x)
Base.:*(n::Integer, a::Ang{T}) where {T} = Ang(n * a.x)
Base.:*(a::Ang{T}, n::Integer) where {T} = Ang(a.x * n)
Base.:-(a::Ang) = Ang(-a.x)

"""
    radians(a::Ang{T})::float(T)

Convert `a` to radians.
"""
radians(a::Ang{T}) where {T} = pi * float(T)(a)

for FT in (:Float64, :Float32)
 @eval function $FT(a::Ang{T}) where {T}
     $FT(a.x / typemax(T))
 end
end

# For now, we use Float64 for everything
function Base.AbstractFloat(a::Ang)
    Float64(a)
end

function Base.:/(a::Ang, n::Integer)
    (q, r) = divrem(a.x, n)
    iszero(r) && return Ang(q)
    x = typeof(q)(q + r/n)
    Ang(x)
end

function Base.show(io::IO, ::PRETTY, a::Ang)
    print(io, Float64(a), " π")
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ang})
    Ang(rand(Int))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ang{T}}) where {T}
    Ang(rand(T))
end

"""
    random_angle(::Type{T}=Float64, shape...)::T where {T}

Return a random angle in `[-pi, pi]`
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

function _random_angle(::Type{IntT}, rng::Random.AbstractRNG, args...) where {IntT}
    radians(rand(rng, Ang{IntT}, args...))
end

function _random_angle(::Type{IntT}, args...) where {IntT}
    radians(rand(Ang{IntT}, args...))
end

function random_phase(::Type{T}=Float64) where {T}
    cis(random_angle(T))
end

for func in (:cos, :sin, :cis, :tan, :sincos, :csc, :sec)
    funcpi = Symbol(func, :pi)
    @eval $func(a::Ang) = $funcpi(Float64(a))
end

end # module Angles
