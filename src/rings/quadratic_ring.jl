module QuadraticRings

import LinearAlgebra
import Base: promote_rule, show, convert
import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two, root_two, coeffs,
    mul_half, conj_root_two, norm_root_two, isrational, isunit, invstrict
import ..Dyadics: Dyadic

import ..Singletons: RootTwoT, RootTwo, Two, InvRootTwo, InvRootTwoT, InvTwo
import ..RootOnes: RootOne
import ..Utils: PRETTY

export QuadraticRing, QuadraticRing2, ZrootD, ZRoot2, DRoot2

########################
####
#### QuadraticRing
####
########################

# QuadraticInteger is a good name if the type of a and by is Integer
# But we allow Dyadic, as well. Not sure what a good
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

If `D` is the type `Dyadic`, then the ring represents `𝔻[√D]`, that is, the ring of
dyadic fractions with `√D` adjoined.

The following aliases are defined:
```julia
const DRoot2{T1, T2} = QuadraticRing{2, Dyadic{T1, T2}}
```

```julia
const ZrootD{T1, T2} = QuadraticRing{T1, T2} where {T1, T2<:Integer}
```

```julia
const ZRoot2{T} = QuadraticRing{2, T} where {T<:Integer}
```

# Examples
```jldoctest
julia> QuadraticRing{2, Dyadic{Int,Int}}(1, 2)
1 + 2√2

julia> repr(QuadraticRing{2, Dyadic{Int,Int}}(1, 2))
"QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(1, 0), Dyadic{Int64, Int64}(2, 0))"
```

# Notes
Values that could be supported, but are not, include integers such that `D < 0` and integers such that `D≡1 (mod 4)`.

The Gaussian integers would be represented by `D` equal to `-1`, if this were supported. But Gaussian
integers are already exactly represented by `Complex{<:Integer}`.
"""
struct QuadraticRing{D, CoeffT <: Real} <: Real
    function QuadraticRing{D, T}(a::T, b::T) where {D, T}
        isa(D, Integer) || throw(ArgumentError(lazy"D must be an Integer"))
        new(a, b)
    end

    a::CoeffT
    b::CoeffT
end

function QuadraticRing{D, CT}(a, b) where {D, CT}
    (a, b) = (CT(a), CT(b))
    return QuadraticRing{D}(a, b)
end

QuadraticRing{D}(a::T, b::T) where {T, D} = QuadraticRing{D, T}(a, b)

"""
    QuadraticRing{D, CT}(x::QuadraticRing) where {D, CT}

Attempt to map an element `x` of one quadratic ring to another.

As pointed out by Ross and Selinger, there is a unique ring homomorphism from `ℤ[√2]` to any ring containing
`√2`. It exists because `ℤ[√2]` is the free such ring. This fact is expressed exactly in their Haskell
implementation. Idiomatic Julia probably does not support this. Or at least it would be an unnecessary complication.
This method for `QuadraticRing` is less restrictive in input and more restrictive in output. Because this is Julia, there is
no guarantee that this construction will succeed for every call that dispatches to this method.
"""
function QuadraticRing{D, CT}(x::QuadraticRing) where {D, CT}
    (a, b) = map(CT, coeffs(x))
    return QuadraticRing{D, CT}(a, b)
end

# Avoid temptation to do D = an integer, QuadraticRing2{D, CoeffT}.
# That kills performance
const QuadraticRing2{CoeffT} = QuadraticRing{2} where CoeffT

"""
    coeffs(q::QuadraticRing)
    coeffs(q::ZRoot2)
    coeffs(q::DRoot2)

Return a `Tuple` of the two coeffients of `q`.

# Examples
```jldoctest
julia> coeffs(ZRoot2(1,2))
(1, 2)

julia> coeffs(DRoot2(1,2))
(Dyadic{Int64, Int64}(1, 0), Dyadic{Int64, Int64}(2, 0))

julia> coeffs(DRoot2(1,Dyadic(3,2)))
(Dyadic{Int64, Int64}(1, 0), Dyadic{Int64, Int64}(3, 2))
```
"""
coeffs(q::QuadraticRing) = (q.a, q.b)

# Can't constrain `D` to be an integer :(
"""
    ZrootD{D, CoeffT} where {D, CoeffT<:Integer}

Represents the ring `ℤ[√D]` of quadratic integers of radix `D`,
provided `D` is a valid integer (See `QuadraticRing`).

`ZrootD` is an alias defined by
```julia
ZrootD{D, CoeffT} = QuadraticRing{D, CoeffT} where {D, CoeffT<:Integer}
```

See also `ZRoot2`, `DRoot2`.

# Examples
```jldoctest
julia> ZrootD{3}(2, 3)
2 + 3√3
```
"""
const ZrootD{T, V} = QuadraticRing{T, V} where {T, V<:Integer}
#const ZrootD{D, T2} = QuadraticRing{D, T2} where {D, T2<:Integer}

function ZrootD{D}(a, b=zero(a)) where D
    (a, b) = promote(a, b)
    ZrootD{D}{typeof(a)}(a, b)
end

"""
    ZRoot2{CoeffT} where {CoeffT<:Integer}

Represents the ring `ℤ[√2]` of quadratic integers of radix `2`.

`ZRoot2` is an alias defined by
```
ZRoot2{CoeffT} = QuadraticRing{2, CoeffT} where {2, CoeffT<:Integer}
```

See also `DRoot2`.

# Examples
```jldoctest
julia> ZRoot2(2, 3)
2 + 3√2

julia> repr(ZRoot2(2, 3))
"QuadraticRing{2, Int64}(2, 3)"

julia> ZRoot2{BigInt}(1, 2)
1 + 2√2

julia> typeof(coeffs(ZRoot2{BigInt}(1, 2)))
Tuple{BigInt, BigInt}

julia> ZRoot2{Dyadic{Int, Int}}(1, 2)  # Enforces ℤ as base ring.
ERROR: TypeError: in QuadraticRing, in T, expected T<:Integer, got Type{Dyadic{Int64, Int64}}
```
"""
const ZRoot2{T} = QuadraticRing{2, T} where {T<:Integer}
# const ZRoot2{T} = ZrootD{2, T}
# Defining in terms of ZrootD may help printing alias type. (Nope)

# This has the effect of printing ZRoot2 as we want.
# But it no longer prints: (alias for ...)
function Base.show(io::IO, ::Type{ZRoot2{T}}) where {T <: Integer}
    print(io, "ZRoot2{$T}")
end

##
## Constructors
##

function ZRoot2(a, b=zero(a))
    (a, b) = promote(a, b)
    ZRoot2{typeof(a)}(a, b)
end


ZRoot2(::RootTwoT) = ZRoot2{Int}(RootTwo)

function QuadraticRing{2, T}(::RootTwoT) where {T}
    QuadraticRing{2, T}(zero(T), one(T))
end

# function ZRoot2{T}(::RootTwoT) where {T}
#     ZRoot2(zero(T), one(T))
# end

isrational(q::QuadraticRing{<:Any, <:Integer}) = iszero(q.b)
isrational(q::QuadraticRing{<:Any, <:Dyadic}) = iszero(q.b)

"""
    root_two(::Type{ZRoot2{T}}) where T

Return a value of `ZRoot2{T}` representing the square root of two.

# Example
```jldoctest
julia> root_two(ZRoot2{Int})
0+ √2
```
"""
root_two(::Type{ZRoot2{T}}) where T = ZRoot2(zero(T), one(T))

function show(io::IO, ::MIME"text/plain", qr::QuadraticRing{D}) where {D}
    if isone(qr.b)
        show(io, MIME"text/plain"(), qr.a)
        print(io, " + √", D)
    else
        show(io, MIME"text/plain"(), qr.a)
        print(io, " + ")
        show(io, MIME"text/plain"(), qr.b)
        print(io, "√", D)
    end
end

# Follow Julia convention by printing both real and imaginary parts
# even if they are == zero.
function show(io::IO, ::MIME"text/plain", qr::Complex{<:QuadraticRing})
    rp = real(qr)
    ip = imag(qr)
    print(io, "(")
    show(io, PRETTY(), rp)
    print(io, ")")
    print(io, " + ")
    print(io, "(")
    show(io, PRETTY(), imag(qr))
    print(io, ")im")
end

"""
    sign(x::QuadraticRing)

Return `cmp(x, zero(x))`
"""
function Base.sign(x::QuadraticRing)
    (a, b) = (x.a, x.b)
    asgn = sign(a)
    bsgn = sign(b)
    iszero(a) && return bsgn
    iszero(b) && return asgn
    (asgn > 0 && bsgn > 0) && return 1
    (asgn < 0 && bsgn < 0) && return -1
    res = cmp(a*a, 2 * b*b)
    asgn == 1 ? res : -res
end

Base.cmp(x::T, y::T) where {T <: QuadraticRing} = sign(x - y)

function Base.:<(x::T, y::T) where {T <: QuadraticRing}
    cmp(x, y) == -1
end

function Base.:<=(x::T, y::T) where {T <: QuadraticRing}
    cmp(x, y) != 1
end

function Base.:(==)(x::T, y::T) where {T <: QuadraticRing}
    x.a == y.a && x.b == y.b
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

function convert(::Type{T}, q::QuadraticRing{D}) where {T<:Number, D}
    (a1, D1, b1) = promote(q.a, D, q.b)
    convert(T, a1) + sqrt(convert(T, D1)) * convert(T, b1)
end

#convert(::Type{T}, ::QuadraticRing{D}) where {T<:Number, D}

function convert(::Type{T}, q::QuadraticRing{D}) where {T<:Integer, D}
    iszero(q.b) || throw(ArgumentError(lazy"Inexact error converting $q to $T"))
    convert(T, q.a)
end

Base.float(q::QuadraticRing) = convert(AbstractFloat, q)

Base.zero(::Type{QuadraticRing{D, CoeffT}}) where {D, CoeffT} = QuadraticRing{D, CoeffT}(zero(CoeffT), zero(CoeffT))
Base.zero(q::QuadraticRing) = zero(typeof(q))
Base.transpose(q::QuadraticRing) = q

function Base.iszero(q::QuadraticRing)
    iszero(q.a) && iszero(q.b)
end

Base.one(::Type{QuadraticRing{D, CoeffT}}) where {D, CoeffT} = QuadraticRing{D, CoeffT}(one(CoeffT), zero(CoeffT))
Base.one(q::QuadraticRing) = one(typeof(q))

"""
    isunit(q::QuadraticRing{<:Any, <:Integer})

Return `true` if `norm_root_two(q)` is either `1` or `-1`.
"""
function isunit(q::QuadraticRing{<:Any, <:Integer})
    nq = norm_root_two(q)
    nq == one(q) || nq == -one(q)
end

invstrict(q::QuadraticRing) = inv(q)

"""
    inv(q::QuadraticRing)

Return `p`, of the same type as `q`, such that `p * q` equals `1`.

If no such `p` exists, throw an error.
"""
function Base.inv(q::QuadraticRing{<:Any, <:Integer})
    nr = norm_root_two(q)
    isone(nr) && return typeof(q)(q.a, -q.b)
    isone(-nr) && return typeof(q)(-q.a, q.b)
    throw(ArgumentError(lazy"Inexact error: $q has no inverse"))
end

function Base.inv(q::QuadraticRing{<:Any, <:Dyadic})
    nr = norm_root_two(q)
    @show nr
    isone(nr) && return typeof(q)(q.a, -q.b)
    isone(-nr) && return typeof(q)(-q.a, q.b)
    throw(ArgumentError(lazy"Inexact error: $q has no inverse"))
end

function Base.big(q::QuadraticRing{D}) where D
    QuadraticRing{2}(big(q.a), big(q.b))
end

function Base.big(q::QuadraticRing{D, T}) where {D, T <: AbstractFloat}
    big(q.a) + sqrt(big(D)) * big(q.b)
end

# This is a misnomer. Works for other D as well
"""
    norm_root_two(qi::QuadraticRing{D})

Return the norm of `qi`.

This method maps `a + b√D` to `a² - D b²`.
This method implements the √2-norm.
"""
norm_root_two(qi::QuadraticRing{D}) where D = qi.a * qi.a  - D * (qi.b * qi.b)

"""
    conj_root_two(qi::QuadraticRing{2})

Maps `a + b√2` to `a - b√2`.

This implements √2-conjugation.
"""
conj_root_two(qi::QuadraticRing{2}) = QuadraticRing{2}(qi.a, -qi.b)

rootDconj(qi::QuadraticRing{D}) where D = QuadraticRing{D}(qi.a, -qi.b)

Base.:*(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D =
    QuadraticRing{D}(q1.a * q2.a + D * q1.b * q2.b, q1.a * q2.b + q1.b * q2. a)

Base.:-(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a - q2.a, q1.b - q2.b)
Base.:-(q::QuadraticRing{D}) where D = QuadraticRing{D}(-q.a, -q.b)
Base.:+(q1::QuadraticRing{D}, q2::QuadraticRing{D}) where D = QuadraticRing{D}(q1.a + q2.a, q1.b + q2.b)
Base.:^(q::QuadraticRing{D}, n::Integer) where D = Base.power_by_squaring(q, n)

function Base.:+(q::QuadraticRing{M, T}, n::Integer) where {M, T}
    (a, b, n) = promote(q.a, q.b, n)
    QuadraticRing{M}(n + a, b)
end

Base.:+(n::Integer, q::QuadraticRing) = q + n

Base.:-(q::QuadraticRing, n::Integer) = typeof(q)(q.a - typeof(q.a)(n), -q.b)
Base.:-(n::Integer, q::QuadraticRing) = typeof(q)(typeof(q.a)(n) - q.a, -q.b)
#Base.:-(n::Integer, q::QuadraticRing) = typeof(q)(typeof(q.a)(n) - q.a, q.b)

Base.:*(q::QuadraticRing{D}, n::Integer) where {D} = QuadraticRing{D}(n * q.a, n * q.b)
Base.:*(n::Integer, q::QuadraticRing{D}) where {D} = q * n

Base.:*(::RootTwoT, q::QuadraticRing{2}) = QuadraticRing{2}(Two * q.b, q.a)
Base.:*(q::QuadraticRing{2}, ::RootTwoT) = RootTwo * q

Base.:*(::RootTwoT, n::Integer) = QuadraticRing{2}(zero(n), n)
Base.:*(n::Integer, ::RootTwoT) = RootTwo * n

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

function canonical(z::Complex{<:QuadraticRing})
    Complex(canonical(real(z)), canonical(imag(z)))
end

promote_rule(::Type{V}, ::Type{T2})  where {V <: QuadraticRing{<:Any, T},T2} where T =
    promote_type(float(T), T2)

"""
    DRoot2{T1, T2}

Represents the ring `𝔻[√2] = ℤ[√½]`.

`T1` is the type of the numerator in the dyadic fractions. `T2` is the type of the exponent in the denominator.

# Examples
```jldoctest
julia> DRoot2(Dyadic(1,2), Dyadic(3, 4))
1/2² + 3/2⁴√2

julia> repr(DRoot2(Dyadic(1,2), Dyadic(3, 4)))
"QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(1, 2), Dyadic{Int64, Int64}(3, 4))"


julia> DRoot2(4, 5)
4 + 5√2

julia> repr(DRoot2(4, 5))
"QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(4, 0), Dyadic{Int64, Int64}(5, 0))"
```
    """
const DRoot2{T, V} = QuadraticRing{2, Dyadic{T, V}}
#const DRoot2{T1, T2} = QuadraticRing{2, Dyadic{T1, T2}}

function DRoot2(a, b)
    a1 = Dyadic(a)
    a2 = Dyadic(b)
    QuadraticRing{2}(a1, a2)
end

"""
    root_two(::Type{DRoot2{T1, T2}}) where {T1, T2}

Return a value of `DRoot2{T1,T2}` representing the square root of two.

# Example
```jldoctest
julia> root_two(DRoot2{Int,Int})
0+ √2
```
"""
function root_two(::Type{DRoot2{T1, T2}}) where {T1, T2}
    z1 = zero(T1)
    z2 = zero(T2)
    o1 = one(T1)
    DRoot2(Dyadic(z1,z2), Dyadic(o1,z2))
end

Base.:*(::InvRootTwoT, q::DRoot2) = QuadraticRing{2}(q.b, InvTwo * q.a)
Base.:*(q::DRoot2, ::InvRootTwoT) = q * InvRootTwo

"""
    one_over_root_two(::Type{DRoot2{T1, T2}}) where {T1, T2}

Return a value representing the reciprocal of the square root of two.

# Example
```jldoctest
julia> one_over_root_two(DRoot2{Int,Int})
0 + 1/2√2

julia> big(one_over_root_two(DRoot2{Int,Int}))
0.0 + 0.5√2

julia> big(big(one_over_root_two(DRoot2{Int,Int})))
0.707106781186547524400844362104849039284835937688474036588339868995366239231051
```
"""
function one_over_root_two(::Type{DRoot2{T1, T2}}) where {T1, T2}
    z1 = zero(T1)
    z2 = zero(T2)
    o1 = one(T1)
    o2 = one(T2)
    DRoot2(Dyadic(z1,z2), Dyadic(o1,o2))
end

function mul_half(q::QuadraticRing{D, <:Dyadic}, n::Integer=1) where {D}
    new_coeffs = map(x -> mul_half(x, n), coeffs(q))
    typeof(q)(new_coeffs...)
end

Base.complex(r::RootOne{8}) = Complex(real(r), imag(r))

function Base.imag(r::RootOne{8})
    k = r.k
    return if k == 0
        DRoot2(0, 0)
    elseif k == 1
        DRoot2(0, Dyadic(1, 1))
    elseif k == 2
        DRoot2(1, Dyadic(0, 0))
    elseif k == 3
        DRoot2(0, Dyadic(1, 1))
    elseif k == 4
        DRoot2(0, 0)
    elseif k == 5
        DRoot2(0, Dyadic(-1, 1))
    elseif k == 6
        DRoot2(-1, Dyadic(0, 0))
    elseif k == 7
        DRoot2(0, Dyadic(-1, 1))
    end
end

function Base.real(r::RootOne{8})
    k = r.k
    return if k == 0
        DRoot2(1, 0)
    elseif k == 1
        DRoot2(0, Dyadic(1, 1))
    elseif k == 2
        DRoot2(0, 0)
    elseif k == 3
        DRoot2(0, Dyadic(-1, 1))
    elseif k == 4
        DRoot2(-1, 0)
    elseif k == 5
        DRoot2(0, Dyadic(-1, 1))
    elseif k == 6
        DRoot2(0, 0)
    elseif k == 7
        DRoot2(0, Dyadic(1, 1))
    end
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(q::QuadraticRing)
        convert($Ti, q)
    end
end

end # module QuadraticRings

#  LocalWords:  Selinger homomorphism
