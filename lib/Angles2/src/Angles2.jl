module Angles2

import Random

export Dar, radtodar, scalepi, unscalepi, minus_one_to_one, radians, random_angle, random_phase

const PRETTY = MIME"text/plain"

"""
    Dar{T}

An angular unit equal to radians scaled by `π`.
"""
struct Dar{T} #  <: Real I think I don't want this
    function Dar(x)
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

Base.:/(d::Dar, x::Real) = Dar(d.x / x)
Base.://(d::Dar, x::Real) = Dar(d.x // x)
Base.div(d::Dar, x::Real) = Dar(div(d.x, x))

Dar(d::Dar) = d

## We might want to use these rotate the angle into a given range.
## For example, some math proofs depend on this.

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

function zero_to_two(x)
    mod(x, 2)
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

# intdiv(x::Integer, n::Integer) = div(x, n)
# intdiv(x::Real, n::Integer) = x / n
# #intdiv(x::Rational, n::Integer) = x * 1//n
# intdiv(x::Dar, n::Integer) = Dar(intdiv(x.x, n))

cscpi(x) = inv(sinpi(x))
secpi(x) = inv(cospi(x))

import Base: cos, sin, cis, tan, sincos, csc, sec

for func in (:cos, :sin, :cis, :tan, :sincos, :csc, :sec)
    funcpi = Symbol(func, :pi)
    @eval $func(dar::Dar) = $funcpi(dar.x)
end

# A random angle between -pi and pi
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

function random_phase(::Type{T}=Float64) where {T}
    cis(random_angle(T))
end

end # module Angles2
