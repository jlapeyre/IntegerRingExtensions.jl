module IntegerExtensions

import LinearAlgebra
import Base: convert
using ILog2

export QuadraticRing, isunit, RootOne, DyadicFraction, CyclotomicRing

########################
####
#### QuadraticRing
####
########################

# QuadraticInteger is a good name if the type of a and by is Integer
# But we allow DyadicFraction, as well. Not sure what a good
# name is. I use QuadraticRing
#
# Assume no k s.t. D != 4k + 1
# i.e. D - 1 not a multiple of 4
# e.g. D == 5 not allowed.
# We do not implement this type for such D
#
# Quadratic integers with D equiv 0 mod 4 do not exist
# We check for neither condition.
struct QuadraticRing{D, CoeffT}
    a::CoeffT
    b::CoeffT
end

QuadraticRing(a::T, b::T, D) where {T} = QuadraticRing{D, T}(a, b)
QuadraticRing{D}(a::T, b::T) where {T, D} = QuadraticRing{D, T}(a, b)

# Maybe not needed
function Base.copy(q::QuadraticRing{D}) where D
    QuadraticRing{D}(q.a, q.b)
end

function convert(::Type{QuadraticRing{D, CoeffT}}, q::QuadraticRing{D}) where {D, CoeffT}
    qa = convert(CoeffT, q.a)
    qb = convert(CoeffT, q.b)
    QuadraticRing{D}(qa, qb)
end

# TODO: Use promote or s.t. like that

convert(::Type{T}, q::QuadraticRing{D}) where {T, D} = convert(T, q.a) + sqrt(D) * convert(T, q.b)

function convert(::Type{T}, q::QuadraticRing{D}) where {T<:Integer, D}
    q.b == 0 || throw(ArgumentError(lazy"Inexact error converting $q to $T"))
    convert(T, q.a)
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(q::QuadraticRing)
        convert($Ti, q)
    end
end

Base.float(q::QuadraticRing) = convert(AbstractFloat, q)
Base.zero(::Type{QuadraticRing{D, CoeffT}}) where {D, CoeffT} = QuadraticRing{D, CoeffT}(0, 0)
Base.zero(q::QuadraticRing) = zero(typeof(q))

"""
    isunit(q::QuadraticRing)

Return `true` if `norm(q)` is either `1` or `-1`.
"""
function isunit(q::QuadraticRing)
    nq = LinearAlgebra.norm(q)
    nq == 1 || nq == -1
end

Base.big(q::QuadraticRing{D}) where D = convert(QuadraticRing{D, BigInt}, q)
LinearAlgebra.norm(qi::QuadraticRing{D}) where D = qi.a * qi.a  - D * (qi.b * qi.b)
Base.conj(qi::QuadraticRing{D}) where D = QuadraticRing{D}(qi.a, -qi.b)

Base.:*(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D =
    QuadraticRing{D}(q1.a * q2.a + 2 * q1.b * q2.b, q1.a * q2.b + q1.b * q2. a)

Base.:-(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a - q2.a, q1.b - q2.b)
Base.:-(q::QuadraticRing{D}) where D = QuadraticRing{D}(-q.a, -q.b)
Base.:+(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a + q2.a, q1.b + q2.b)
Base.:^(q::QuadraticRing{D}, n::Integer) where D = Base.power_by_squaring(q, n)

########################
####
#### RootOne
####
########################

# This struct implements the cyclic group of order N,
# plus a few methods particular to the interpretation as roots of unity.
"""
    RootOne{N}

`N`th roots of unity
"""
struct RootOne{N}
    function RootOne{N}(k) where {N}
        new{N}(mod(k, N))
    end
    k::Int # Could use smaller `k`.
end

# Examples:
# convert(ComplexF64, z)
# convert(Complex{BigFloat}, z)
function convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N}
    kT = convert(T, r.k)
    cispi(2 * kT / N)
end

RootOne(k, N) = RootOne{N}(k)
Base.one(::Type{RootOne{N}}) where {N} = RootOne{N}(0)
Base.one(x::RootOne) = one(typeof(x))
Base.float(r::RootOne) = complex(r)
Base.complex(r::RootOne) = convert(ComplexF64, r)
Base.big(r::RootOne) = convert(Complex{BigFloat}, r)
# TODO: angle should agree with floating point angle: in [-pi, pi] rather than [0, 2pi]
Base.angle(r::RootOne{N}) where {N} = 2 * pi * r.k / N
Base.abs2(r::RootOne) = 1
Base.abs(r::RootOne) = 1
Base.real(r::RootOne{N}) where {N} = cospi(2 * r.k / N)
Base.imag(r::RootOne{N}) where {N} = sinpi(2 * r.k / N)
Base.isreal(r::RootOne{N}) where {N} = r.k == 0 || div(N, r.k) == 2
Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k + r2.k)
Base.:^(r::RootOne{N}, n::Integer) where {N} = RootOne{N}(r.k * n)
Base.inv(r::RootOne{N}) where {N} = RootOne{N}(N - r.k)

function Base.:-(r::RootOne{N}) where {N}
    iseven(N) || throw(ArgumentError(lazy"InexactError unary minus of type RootOne{$N}"))
    RootOne{N}(r.k + div(N, 2))
end

########################
####
#### DyadicFraction
####
########################

struct DyadicFraction{aT, kT}
    a::aT
    k::kT
end

function convert(::Type{T}, f::DyadicFraction) where {T}
    if iszero(f.k)
        return convert(T, f.a)
    end
    convert(T, f.a) / convert(T, 2)^f.k
end

function convert(::Type{Float64}, f::DyadicFraction)
    if iszero(f.k)
        return convert(Float64, f.a)
    end
    convert(Float64, f.a) * 0.5^f.k
end

function convert(::Type{T}, f::DyadicFraction) where {T <: Integer}
    f.k == 0 || throw(ArgumentError(lazy"Inexact error converting $f to $T"))
    convert(T, f.a)
end

function convert(::Type{DyadicFraction}, n::Integer)
    DyadicFraction(n, zero(n))
end

function convert(::Type{DyadicFraction{T, K}}, n::T) where {T <: Integer, K}
    DyadicFraction{T, K}(n, 0)
end

convert(::Type{Rational{T}}, f::DyadicFraction) where {T} =
    Rational{T}(convert(T, f.a), convert(T, 2)^f.k)

DyadicFraction(r::Rational) = convert(DyadicFraction, r)
DyadicFraction(n::Integer) = DyadicFraction(n, zero(n))

Base.Rational(f::DyadicFraction{aT}) where {aT} = convert(Rational{aT}, f)

function convert(::Type{DyadicFraction}, r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    DyadicFraction(r.num, ilog2(r.den))
end

Base.float(f::DyadicFraction) = convert(Float64, f)
Base.big(f::DyadicFraction) = convert(BigFloat, f)
Base.:*(f1::DyadicFraction, f2::DyadicFraction) = DyadicFraction(f1.a * f2.a, f1.k + f2.k)
Base.:-(f::DyadicFraction) = DyadicFraction(-f.a, f.k)
Base.:*(n::Integer, f::DyadicFraction) = DyadicFraction(n * f.a, f.k)

Base.:+(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, +)
Base.:-(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, -)
function _plus(f1::DyadicFraction, f2::DyadicFraction, op)
    (minex, maxex) = minmax(f1.k, f2.k)
    DyadicFraction(op(2^(f2.k - minex) * f1.a,  2^(f1.k - minex) * f2.a), maxex)
end

########################
####
#### CyclotomicRing
####
########################

struct CyclotomicRing{M, CoeffT}
    coeffs::NTuple{M, CoeffT}
end

function Base.:+(c1::CyclotomicRing{M, CoeffT}, c2::CyclotomicRing{M, CoeffT}) where {M, CoeffT}
    CyclotomicRing{M, CoeffT}(c1.coeffs .+ c2.coeffs)
end

end # module IntegerExtensions
