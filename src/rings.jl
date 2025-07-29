module Rings

import LinearAlgebra
import Base: convert, zero, iszero, one, isone
import ..IntegerExtensions: imaginary, sqrt_imaginary, one_over_root_two, canonical
using ILog2

export QuadraticRing, isunit, RootOne, DyadicFraction, CyclotomicRing

export Domega

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

# Avoid temptation to do D = an integer, QuadraticRing2{D, CoeffT}.
# That kills performance
const QuadraticRing2{CoeffT} = QuadraticRing{2}

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
    iszero(q.b) || throw(ArgumentError(lazy"Inexact error converting $q to $T"))
    convert(T, q.a)
end


Base.float(q::QuadraticRing) = convert(AbstractFloat, q)
Base.zero(::Type{QuadraticRing{D, CoeffT}}) where {D, CoeffT} = QuadraticRing{D, CoeffT}(zero(CoeffT), zero(CoeffT))
Base.zero(q::QuadraticRing) = zero(typeof(q))

function Base.iszero(q::QuadraticRing)
    iszero(q.a) && iszero(q.b)
end

Base.one(::Type{QuadraticRing{D, CoeffT}}) where {D, CoeffT} = QuadraticRing{D, CoeffT}(one(CoeffT), zero(CoeffT))
Base.one(q::QuadraticRing) = one(typeof(q))

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
    struct RootOne{N}
    RootOne{N}(k::Integer)

`N`th roots of unity

`N` should always be a literal or a `const`. Otherwise performance
is severely degraded.

The type `const RootOne8 = RootOne{8}` is defined.
`k` will be stored as `mod(k, N)`, which takes values in `(0,...,N-1)`.
# Examples
```jldoctest

julia> RootOne8(1)
RootOne8(1)
```

```jldoctest
julia> RootOne8(-1)
RootOne8(7)

julia> -RootOne8(1) # Unary minus
RootOne8(5)

julia> -RootOne8(5)
RootOne8(1)

julia> -RootOne{9}(1)
ERROR: ArgumentError: InexactError unary minus of type RootOne{9}

julia> RootOne8(1)^8
RootOne8(0)

julia> isone(RootOne8(1)^8)
true

julia> one(RootOne8)
RootOne8(0)

julia> one(RootOne8(5))
RootOne8(0)

julia> RootOne8(2) * RootOne8(3)
RootOne8(5)
```
"""
struct RootOne{N}
    function RootOne{N}(k) where {N}
        new{N}(mod(k, N))
    end
    k::Int # Could use smaller `k`.
end

# This is largely for performance. If N is not constant,
# performance suffers by factor of 100X or more.
const RootOne8 = RootOne{8}

"""
    imaginary(::Type{RootOne8})

The value of type `RootOne8` that represents the imaginary unit.

This is `RootOne8(2)`.
"""
function imaginary(::Type{RootOne8})
    RootOne8(2)
end

"""
    sqrt_imaginary(::Type{RootOne8})

The value of type `RootOne8` that represents the principal square root of the imaginary unit.

This is `RootOne8(1)`.
"""
function sqrt_imaginary(::Type{RootOne8})
    RootOne8(1)
end

# If N is not a compile-time constant in RootOne{N}, then
# operations are very, very, inefficient.
# Either make N constant, eg. `const N = 8`, or
# `const RootOne8 = RootOne{8}`
# Otherwise use a type like RootOneA with no type parameter.
#
# Same as RootOne{N}, but don't use type param for N
struct RootOneA
    function RootOneA(k, N)
        new(mod(k, N), N)
    end
    k::Int
    N::Int
end

Base.:*(r1::RootOneA, r2::RootOneA) = RootOneA(r1.k + r2.k, r1.N)

# Examples:
# convert(ComplexF64, z)
# convert(Complex{BigFloat}, z)
function convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N}
    Complex{T}(r)
end

function Complex{T}(r::RootOne{N}) where {T, N}
    kT = convert(T, r.k)
    cispi(2 * kT / N)
end

"""
    RootOne(k, N)

Construct `RootOne{N}(k)`.
"""
RootOne(k, N) = RootOne{N}(k)
Base.one(::Type{RootOne{N}}) where {N} = RootOne{N}(0)
Base.one(x::RootOne) = one(typeof(x))

Base.float(r::RootOne) = complex(r)
Base.complex(r::RootOne) = convert(ComplexF64, r)
Base.big(r::RootOne) = convert(Complex{BigFloat}, r)
# TODO: angle should agree with floating point angle: in [-pi, pi] rather than [0, 2pi]

function Base.angle(r::RootOne{N}) where {N}
    k = r.k <= N/2 ? r.k : r.k - N
    2 * pi * k / N
end

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

"""
    canonical(df::DyadicFraction)

Return `df` in canonical form. That is, with `df.k` as small as possible.
"""
function canonical(df::DyadicFraction)
    iszero(df.a) && return DyadicFraction(df.a, zero(df.k))
    (c, m) = factortwos(df.a)
    iszero(c) && return df
    k = df.k
    if c > k
        return typeof(df)(2^(c-k)*m, 0)
    elseif c < k
        return typeof(df)(m, k - c)
    else
        return typeof(df)(m, 0)
    end
end

"""
    factortwos(n)

Factor `n` as `m * 2^c` and return `(c, m)`.
"""
function factortwos(n)
    n == 0 && return (0, 0)
    c = 0
    m = n
    while true
        isodd(m) && return (c, m)
        m = div(m, 2)
        c += 1
    end
end

function zero(::Type{DyadicFraction{aT, bT}}) where {aT, bT}
    DyadicFraction(zero(aT), zero(bT))
end

# Careful. There is more than on way to represent zero.
function iszero(df::DyadicFraction)
    iszero(df.a)
end

function one(::Type{DyadicFraction{aT, bT}}) where {aT, bT}
    DyadicFraction(one(aT), zero(bT))
end

one(::DyadicFraction{aT, bT}) where {aT, bT} = one(DyadicFraction{aT, bT})

function canonical(qr::QuadraticRing)
    typeof(qr)(canonical(qr.a), canonical(qr.b))
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
    iszero(f.k) && return convert(T, f.a)
    cf = canonical(f)
    cf.k == 0 || throw(ArgumentError(lazy"Inexact error converting $f to $T"))
    convert(T, cf.a)
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

function Complex{Tc}(df::DyadicFraction) where {Tc}
    T = Complex{Tc}
    T(df.a) * T(1//2)^df.k
end

function Base.:(==)(df1::DyadicFraction, df2::DyadicFraction)
    if df1.k > df2.k
        return df1.a == 2^(df1.k - df2.k) * df2.a
    elseif df1.k < df2.k
        return df2.a == 2^(df2.k - df1.k) * df1.a
    else
        return df1.a == df2.a
    end
end

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
Base.:*(f::DyadicFraction, n::Integer) = n * f

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

# The type of the exponent of (1/2) is hardcoded to Int.
# We probably only need one type for this field.
const Domega{T} = CyclotomicRing{4, DyadicFraction{T, Int}}

"""
    imaginary(::Type{Domega{T}}) where {T}

The value of type `Domega{T}` that represents the imaginary unit.
"""
function imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
    CyclotomicRing(z, o, z, z)
end

"""
    sqrt_imaginary(::Type{Domega{T}}) where {T}

The value of type `Domega{T}` that represents the principal square root of the imaginary unit.
    """
function sqrt_imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
    CyclotomicRing(z, z, o, z)
end

function one_over_root_two(::Type{CyclotomicRing{4, DyadicFraction{T, Int}}}) where {T}
    DT = DyadicFraction{T, Int}
    z = zero(DT)
    half = DyadicFraction(T(1), 1)
    CyclotomicRing(-half, z, half, z)
end

CyclotomicRing(coeffs::T...) where T = CyclotomicRing(coeffs)

function canonical(c::CyclotomicRing)
    CyclotomicRing(map(canonical,  c.coeffs))
end

function Base.:+(c1::CyclotomicRing{M, CoeffT}, c2::CyclotomicRing{M, CoeffT}) where {M, CoeffT}
    CyclotomicRing{M, CoeffT}(c1.coeffs .+ c2.coeffs)
end

function Base.:-(c1::CyclotomicRing{M, CoeffT}, c2::CyclotomicRing{M, CoeffT}) where {M, CoeffT}
    CyclotomicRing{M, CoeffT}(c1.coeffs .- c2.coeffs)
end

function Base.:-(c::CyclotomicRing)
    CyclotomicRing(.- c.coeffs)
end

function Base.:*(c1::CyclotomicRing{4, CoeffT}, c2::CyclotomicRing{4, CoeffT}) where {CoeffT}
    (a1, b1, c1, d1) = c1.coeffs
    (a2, b2, c2, d2) = c2.coeffs
    coeffs = (
        a2*d1 + c1*b2 + b1*c2 + a1*d2,
        b2*d1 + c1*c2 + b1*d2 - a1*a2,
        c2*d1 + c1*d2 - b1*a2 - a1*b2,
        d1*d2 - b1*b2 - a1*c2 - a2*c1
    )
    CyclotomicRing(coeffs)
end

function Base.one(::Type{CyclotomicRing{4, CoeffT}}) where {CoeffT}
    z = zero(CoeffT)
    o = one(CoeffT)
    CyclotomicRing(z, z, z, o)
end

Base.one(::CyclotomicRing{4, CoeffT}) where {CoeffT} = one(CyclotomicRing{4, CoeffT})

function Base.zero(::Type{CyclotomicRing{4, CoeffT}}) where {CoeffT}
    z = zero(CoeffT)
    CyclotomicRing(z, z, z, z)
end

Base.zero(cr::CyclotomicRing) = zero(typeof(cr))

function Base.convert(::Type{CyclotomicRing{M, CT1}}, c::CyclotomicRing{M, CT2}) where {M, CT1, CT2}
    CyclotomicRing{M, CT1}(c)
end

function CyclotomicRing{M, CT1}(c::CyclotomicRing{M, CT2}) where {M, CT1, CT2}
    coeffs = convert.(CT1, c.coeffs)
    CyclotomicRing(coeffs)
end

function Base.AbstractFloat(c::CyclotomicRing{4})
    (a, b, c, d) = c.coeffs
    float(d) + float(c) * float(RootOne8(1)) + float(b) * float(RootOne8(2)) +
        float(a) * float(RootOne8(3))
end

function Base.Complex{Tc}(c::CyclotomicRing{4}) where Tc
    (a, b, c, d) = c.coeffs
    T = Complex{Tc}
    T(d) + T(c) * T(RootOne8(1)) + T(b) * T(RootOne8(2)) +
        T(a) * T(RootOne8(3))
end

function Base.big(c::CyclotomicRing)
    Complex{BigFloat}(c)
end

function Base.convert(::Type{Complex{Tc}}, c::CyclotomicRing) where Tc
    Complex{Tc}(c)
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(q::QuadraticRing)
        convert($Ti, q)
    end
    @eval function (::Type{Base.$Ti})(df::DyadicFraction)
        convert($Ti, df)
    end
end

end # module Rings
