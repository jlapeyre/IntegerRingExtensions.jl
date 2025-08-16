module Angles

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

Base.AbstractFloat(::Type{Dar{T}}) where {T} = float(T)
Base.AbstractFloat(a::Dar{T}) where {T} = unscalepi(float(a.x))
Base.Float64(a::Dar) = unscalepi(float(a.x))

Base.convert(::Type{Float64}, a::Dar) = Float64(a)

Base.one(::Type{Dar{T}}) where {T} = one(T)
Base.one(::Dar{T}) where {T} = one(T)

Base.zero(::Type{Dar{T}}) where {T} = Dar(zero(T))
Base.zero(::Dar{T}) where {T} = Dar(zero(T))

Base.:+(d::Dar, x::Real) = Dar(d.x + x)
Base.:-(d::Dar, x::Real) = Dar(d.x - x)
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

# # The following functions work. But they are far
# # From optimal. I stopped when I finally go them
# # to pass some tests.
# function zero_to_two(x)
#     if x >= 0  # There must be a Base function call for this.
#         (q, r) = divrem(x, 2)
#         res =
#             if iszero(q)
#                 r
#             elseif iszero(r)
#                 2 * one(r)
#             else
#                 r
#             end
#         (0 <= res <= 2) || error(lazy"sad outcome $ret")
#         return res
#     else
#         n = round(x)
#         y = x - 2 * n + 2
#         y > 0 || error(lazy"I dont know what to do with $x")
#         zero_to_two(y)
#     end
# end

# function _minus_one_to_one(x)
#     y = zero_to_two(x)
#     if 0 <= y <= 1
#         return y
#     end
#     return y - 2
# end

# function old_minus_one_to_one(x)
#     y = _minus_one_to_one(x)
#     (-1 <= y <= 1) || error(lazy"Sad case of $y")
#     return y
# end

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

function dartorad(dar::Dar)
    dar.x * pi
end

scalepi(x::Real) = x / pi
scalepi(x::Dar) = x
unscalepi(x::Real) = x * pi
unscalepi(x::Dar) = x

intdiv(x::Real, n::Integer) = x / n
intdiv(x::Rational, n::Integer) = x * 1//n
intdiv(x::Dar, n::Integer) = Dar(intdiv(x.x, n))

# function dartorad(x::Real)
#     dar.x * pi
# end

import Base: cos, sin, cis, tan, sincos  # , cospi, sinpi, cispi, tanpi

for func in (:cos, :sin, :cis, :tan, :sincos)
    funcpi = Symbol(func, :pi)
    @eval $func(dar::Dar) = $funcpi(dar.x)
end

Base.isapprox(a::Dar, b::Dar; kws...) = isapprox(a.x, b.x; kws...)
Base.isapprox(a::Dar, x::Real; kws...) = isapprox(unscalepi(a.x), x; kws...)
Base.isapprox(x::Real, a::Dar ; kws...) = isapprox(a, x; kws...)

struct Ang{T<:Integer}
    x::T
end

function Ang(z::Float64)
    zp = minus_one_to_one(z)
    if zp == 1
        Ang(typemax(Int))
    elseif zp == -1
        Ang(-typemax(Int))
    else
        Ang(Int(zp * typemax(Int)))
    end
end

# function Ang{T}(z::Float64) where {T<:Integer}
#     Ang(T(minus_one_to_one(z) * typemax(T)))
# end

function Base.:(==)(a::Ang{T},  b::Ang{T}) where {T}
    a.x == b.x
end

Base.:+(a::Ang{T}, b::Ang{T}) where {T} = Ang(a.x + b.x)
Base.:-(a::Ang{T}, b::Ang{T}) where {T} = Ang(a.x - b.x)
Base.:*(n::Integer, a::Ang{T}) where {T} = Ang(n * a.x)
Base.:*(a::Ang{T}, n::Integer) where {T} = Ang(a.x * n)

for FT in (:Float64, :Float32)
 @eval function $FT(a::Ang{T}) where {T}
     $FT(a.x / typemax(T))
 end
end

function Base.show(io::IO, ::PRETTY, a::Ang)
    print(io, Float64(a), " π")
end

for func in (:cos, :sin, :cis, :tan, :sincos)
    funcpi = Symbol(func, :pi)
    @eval $func(a::Ang) = $funcpi(Float64(a))
end


end # module Angles

