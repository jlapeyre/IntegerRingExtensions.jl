module QuadraticRings

import LinearAlgebra
import Base: promote_rule, show, convert
import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two
import ..DyadicFractions: DyadicFraction

export QuadraticRing, QuadraticRing2, ZrootD, Zroot2, Droot2

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

See also `Zroot2`, `Droot2`.
"""
const ZrootD{T1, T2} = QuadraticRing{T1, T2} where {T1, T2<:Integer}

"""
    Zroot2{CoeffT} where {CoeffT<:Integer}

Represents the ring `ℤ[√2]`

`Zroot2` is an alias defined by
```
Zroot2{CoeffT} = QuadraticRing{2, CoeffT} where {2, CoeffT<:Integer}
```

See also `Droot2`.
"""
const Zroot2{T} = QuadraticRing{2, T} where {T<:Integer}

QuadraticRing(a::T, b::T, D) where {T} = QuadraticRing{D, T}(a, b)
QuadraticRing{D}(a::T, b::T) where {T, D} = QuadraticRing{D, T}(a, b)

function show(io::IO, ::MIME"text/plain", qr::QuadraticRing{D}) where {D}
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

function canonical(qr::QuadraticRing)
    typeof(qr)(canonical(qr.a), canonical(qr.b))
end

promote_rule(::Type{V}, ::Type{T2})  where {V <: QuadraticRing{<:Any, T},T2} where T =
    promote_type(float(T), T2)

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

function one_over_root_two(::Type{Droot2{T1, T2}}) where {T1, T2}
    z1 = zero(T1)
    z2 = zero(T2)
    o1 = one(T1)
    o2 = one(T2)
    Droot2(DyadicFraction(z1,z2), DyadicFraction(o1,o2))
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(q::QuadraticRing)
        convert($Ti, q)
    end
end


end # module QuadraticRings
