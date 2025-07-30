module Rings

import LinearAlgebra
import Base: convert, zero, iszero, one, isone, promote_rule
import ..IntegerExtensions: imaginary, sqrt_imaginary, one_over_root_two, canonical
import ..IntegerExtensions.Utils: subscript, superscript
import ..IntegerExtensions.RootOnes: RootOne8


using ILog2

export QuadraticRing, isunit, DyadicFraction, CyclotomicRing, rootDconj

export Domega, Droot2, ZrootD

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
"""
    QuadraticRing{D, CoeffT}

Represents a ring formed by adjoining the square root of the integer `D` to the ring represented by `CoeffT`.
Integer values of `D` that are equal to 0 or 1 (mod 4) will give incorrect results. Examples of good
values of `D` are `2` and `3`.

If `D` is a supported integer and `CoeffT <: Integer`, then the type represents `ℤ[√D]`, a ring of quadratic integers.

If `D` is the type `DyadicFraction`, then the ring represents `𝔻[√D]`, that is, the ring if dyadic fractions with
`√D` adjoined. For example, the alias `Droot2` is defined as
```julia
Droot2{T1, T2} = QuadraticRing{2, DyadicFraction{T1, T2}}
```

# Notes
Values that could be supported, but are not, include integers such that `D < 0` and integers such that `D≡1 (mod 4)`.

The Gaussian integers would be represented by `D` equal to `-1`, if this were supported. But Gaussian
integers are already exactly represented by `Complex{<:Integer}`.
"""
struct QuadraticRing{D, CoeffT}
    a::CoeffT
    b::CoeffT
end

# Avoid temptation to do D = an integer, QuadraticRing2{D, CoeffT}.
# That kills performance
const QuadraticRing2{CoeffT} = QuadraticRing{2}

# Can't constrain `D` to be an integer :(
"""
    ZrootD{D, CoeffT} where {D, CoeffT<:Integer}

Represents the ring `ℤ[√D]`, provided `D` is a valid integer (See `QuadraticRing`).

`ZrootD` is an alias defined by
```
ZrootD{D, CoeffT} = QuadraticRing{D, CoeffT} where {D, CoeffT<:Integer}
```

See also `Droot2`.
"""
const ZrootD{D, CoeffT} = QuadraticRing{D, CoeffT} where {D, CoeffT<:Integer}

QuadraticRing(a::T, b::T, D) where {T} = QuadraticRing{D, T}(a, b)
QuadraticRing{D}(a::T, b::T) where {T, D} = QuadraticRing{D, T}(a, b)

function Base.show(io::IO, ::MIME"text/plain", qr::QuadraticRing{D}) where {D}
    if isone(qr.b)
        show(io, MIME"text/plain"(), qr.a)
        print(io, "+ √", D)
    else
        show(io, MIME"text/plain"(), qr.a)
        print(io, " + ")
        show(io, MIME"text/plain"(), qr.b)
        print(io, "√", D)
#        print(io, qr.a, " + ", qr.b, "√", D)
    end
end

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

function Base.big(q::QuadraticRing{D}) where D
    QuadraticRing{2}(big(q.a), big(q.b))
end

function Base.big(q::QuadraticRing{D, T}) where {D, T <: AbstractFloat}
    big(q.a) + sqrt(big(D)) * big(q.b)
end

LinearAlgebra.norm(qi::QuadraticRing{D}) where D = qi.a * qi.a  - D * (qi.b * qi.b)

root2conj(qi::QuadraticRing{2}) = QuadraticRing{2}(qi.a, -qi.b)

rootDconj(qi::QuadraticRing{D}) where D = QuadraticRing{D}(qi.a, -qi.b)

Base.:*(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D =
    QuadraticRing{D}(q1.a * q2.a + D * q1.b * q2.b, q1.a * q2.b + q1.b * q2. a)

Base.:-(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a - q2.a, q1.b - q2.b)
Base.:-(q::QuadraticRing{D}) where D = QuadraticRing{D}(-q.a, -q.b)
Base.:+(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a + q2.a, q1.b + q2.b)
Base.:^(q::QuadraticRing{D}, n::Integer) where D = Base.power_by_squaring(q, n)

function Base.:(==)(q::QuadraticRing, n::Integer)
    return iszero(q.b) && n == q.a
end

function Base.:(==)(q::QuadraticRing, x::AbstractFloat)
    (q1, x1) = promote(q, x)
    q1 == x1
end

promote_rule(::Type{V}, ::Type{T2})  where {V <: QuadraticRing{<:Any, T},T2} where T =
    promote_type(float(T), T2)

########################
####
#### DyadicFraction
####
########################

"""
    DyadicFraction{aT, kT}

Represents the ring `𝔻 = ℤ[½]`.

`aT` is the type of the numerator.
`bT` is the type of the exponent on `2` in the denominator.

# Examples
```jldoctest
julia> DyadicFraction(5, 3)
5/2³
```
"""
struct DyadicFraction{aT, kT}
    a::aT
    k::kT
end

"""
    Droot2{T1, T2}

Represents the ring `𝔻[√2] = ℤ[1/√2]`.

`T1` is the type of the numerator in the dyadic fractions. `T2` is the type of the exponent in the denominator.

# Examples
```jldoctest
julia> Droot2(DyadicFraction(1,2), DyadicFraction(3, 4))
1/2² + 3/2⁴√2

julia> Droot2(4, 5) # These integers are still stored as `DyadicFraction`s
4 + 5√2
```
"""
const Droot2{T1, T2} = QuadraticRing{2, DyadicFraction{T1, T2}}

function Droot2(a, b)
    a1 = DyadicFraction(a)
    a2 = DyadicFraction(b)
    QuadraticRing{2}(a1, a2)
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{DyadicFraction{T2,V2}}) where {T1,T2,V1,V2}
    T = promote_type(T1, T2)
    V = promote_type(V1, V2)
    DyadicFraction{T, V}
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:Integer}
    T3 = promote_type(T1, T)
    DyadicFraction{T3, V1}
end

function Base.:^(df::DyadicFraction, n::Integer)
    # Can't copy Int, but really should copy BigInt
    # Need logic for this??
    n == 1 && return df
    return Base.power_by_squaring(df, n)
end

function promote_rule(::Type{DyadicFraction{Int,Int}}, ::Type{Int})
    DyadicFraction{Int, Int}
end

function Base.show(io::IO, ::MIME"text/plain", df::DyadicFraction)
#    print(io, df.a, " / 2^", df.k)
    #    print(io, df.a, " / 2", superscript(df.k))
    if iszero(df.a)
        print(io, zero(df.a))
    else
        if iszero(df.k)
            print(io, df.a)
        else
            print(io, df.a, "/2")
            if !isone(df.k)
               print(io, superscript(df.k))
            end
        end
    end
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

zero(::DyadicFraction{aT, bT}) where {aT, bT} = zero(DyadicFraction{aT, bT})

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

function convert(::Type{DyadicFraction{T,V}}, f::DyadicFraction) where {T, V}
    DyadicFraction{T, V}(T(f.a), V(f.k))
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

convert(::Type{DyadicFraction{T, K}}, n::T1) where {T <: Integer, T1 <: Integer, K} = DyadicFraction{T, K}(n)

function DyadicFraction{T, K}(n::T1)  where {T <: Integer, T1 <: Integer, K}
    DyadicFraction{T, K}(convert(T, n), 0)
end

convert(::Type{Rational{T}}, f::DyadicFraction) where {T} =
    Rational{T}(convert(T, f.a), convert(T, 2)^f.k)

DyadicFraction(r::Rational) = convert(DyadicFraction, r)
DyadicFraction(n::Integer) = DyadicFraction(n, zero(n))
DyadicFraction(x::DyadicFraction) = x
DyadicFraction{T,V}(x::DyadicFraction{T,V}) where {T,V} = x

Base.Rational(f::DyadicFraction{aT}) where {aT} = convert(Rational{aT}, f)

function convert(::Type{DyadicFraction}, r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    DyadicFraction(r.num, ilog2(r.den))
end

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

function Base.show(io::IO, ::MIME"text/plain", cr::CyclotomicRing)
    c = cr.coeffs
    n = length(c)
    showcount = 0
    for i in 1:n
        iszero(c[i]) && continue
        showcount += 1
        if showcount > 1
            print(io, " + ")
        end
        if isone(-c[i])
            print(io, "-")
        elseif !isone(c[i])
            show(io, MIME"text/plain"(), c[i])
        end
        if isone(n - i)
            print(io, " ω")
        else
            print(io, " ω", superscript(n-i))
        end
    end
    if showcount == 0
        print(io, 0)
        # TODO, use zero of first(c) somehow.
        # This will be more robust
#        show(io, MIME"text/plain", zero(first(c)))
    end
end

# The type of the exponent of (1/2) is hardcoded to Int.
# We probably only need one type for this field.
"""
    Domega{T}

Represents the ring `𝔻[ω] = ℤ[1/√2, i]`.

The type `T<:Integer` is the type of the numerator in the dyadic fractions.

`Domega{T}` is defined as the alias
```
Domega{T} = CyclotomicRing{4, DyadicFraction{T, Int}}
```
"""
const Domega{T} = CyclotomicRing{4, DyadicFraction{T, Int}}

function Base.conj(cyc::Domega{T}) where T
    (a, b, c, d) = cyc.coeffs
    # I think promotion not needed.
    newcoeff = promote(-c, -b, -a, d)
    Domega{T}(newcoeff)
end

function Base.conj(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    newcoeff = (-c, -b, -a, d)
    typeof(cyc)(newcoeff)
end

Base.adjoint(cyc::CyclotomicRing) = conj(cyc)

function root2conj(cyc::Domega{T}) where T
    (a, b, c, d) = cyc.coeffs
    newcoeff = (-a, b, -c, d)
    Domega{T}(newcoeff)
end

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

CyclotomicRing(coeffs::T...) where T = CyclotomicRing{length(coeffs), T}(coeffs)

function CyclotomicRing(c1, coeffs...)
    cs = promote(c1, coeffs...)
    CyclotomicRing(cs)
end

function CyclotomicRing{4, T}(a, b, c, d) where T
    cs = (T(a),T(b),T(c),T(d))
    CyclotomicRing(cs)
end

# TODO: Promote
function Domega(a, b, c, d)
    coeffs = promote(a, b, c, d)
    D = DyadicFraction
    coeffs1 = map(D, coeffs)
    T = typeof(coeffs1[1])
    CyclotomicRing{4, T}(coeffs...)
end

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

function Base.:*(c::CyclotomicRing, r::Real)
    cs = c.coeffs .* r
    CyclotomicRing(cs)
end

Base.:*(r::Real, c::CyclotomicRing) = c * r

Base.:(==)(c1::CyclotomicRing, c2::CyclotomicRing) = c1.coeffs == c2.coeffs

function Base.:^(c::CyclotomicRing, n::Integer)
    n == 0 && return one(c)
    n == 1 && return c
    n == 2 && return c * c
    return Base.power_by_squaring(c, n)
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
