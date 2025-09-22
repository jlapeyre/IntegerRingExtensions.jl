# @stable module CyclotomicRings
module CyclotomicRings

#using DispatchDoctor: @stable, @unstable

import LinearAlgebra
import Base: convert, zero, one, promote_rule
import IsApprox: AbstractApprox, Equal

import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two, root_two, coeffs,
    mul_root_two, mul_one_over_root_two, mul_half, conj_root_two, mul_two
import ..Utils: superscript, iszero_strong, isone_strong, PRETTY
import ..RootOnes: RootOne
import ..Dyadics: Dyadic
import ..QuadraticRings: DRoot2, ZRoot2
import ..Singletons: InvTwo, InvTwoT,
    RootTwo, RootTwoT, Two, TwoT,
    InvRootTwo, InvRootTwoT,
    Imag, ImagT,
    RootImag, RootImagT,
    Pow

export CyclotomicRing, ZOmega, DOmega

########################
####
#### CyclotomicRing
####
########################

"""
    struct CyclotomicRing{M, CoeffT}
        coeffs::NTuple{M, CoeffT}
    end

A super-ring of the ring of cyclotomic integers of degree `2M`.

If `CoeffT <: Integer`, then this is represents the ring of cyclotomic integers `ℤ[ω]`.

The symbol `ω` is used when displaying values. It is understood to be the prinicpal `2M`th root of unity.
(See `RootOne{D}`.)

Aliases for particular rings are:

The ring of cyclotomic integers `ℤ[ω]`:
```julia
const ZOmega{T} = CyclotomicRing{4, T} where {T <: Integer}
```

The ring `𝔻[ω] = ℤ[1/√2, i]`:
```julia
const DOmega{T} = CyclotomicRing{4, Dyadic{T, Int}}
```
"""
struct CyclotomicRing{M, CoeffT <: Number} <: Number
    function CyclotomicRing{M, CoeffT}(coeffs::NTuple{M, CoeffT}) where {M, CoeffT <: Number}
        new{M, CoeffT}(coeffs)
    end
    function CyclotomicRing{M, CoeffT}(coeffs::NTuple{M, CoeffT}) where {M, CoeffT <: Union{AbstractFloat, Complex{<:AbstractFloat}}}
        throw(ArgumentError(lazy"Coefficients of type $CoeffT not supported."))
    end
    coeffs::NTuple{M, CoeffT}
end

# Thank you Elrod
# https://discourse.julialang.org/t/is-this-test-detect-unbound-args-result-valid-or-a-bug/96987
function CyclotomicRing(coeffs::Tuple{T,Vararg{T,nMinus1}}) where {T<:Number, nMinus1}
    CyclotomicRing{nMinus1 + 1, T}(coeffs)
end

###
### Aliases for special cases: DOmega, ZOmega
###

"""
    DOmega{T}

Represents the ring `𝔻[ω] = ℤ[1/√2, i]`.

The type `T<:Integer` is the type of the numerator in the dyadic fractions.

Here, `𝔻 = ℤ[½]` is the ring of dyadic fractions, implemented by `Dyadic`.

`DOmega{T}` is defined as the alias

```julia
DOmega{T} = CyclotomicRing{4, Dyadic{T, Int}}
```
"""
const DOmega{T} = CyclotomicRing{4, Dyadic{T, Int}}

"""
    ZOmega{T <: Integer}

Represents the ring `ℤ[ω]`.

The type `T<:Integer` is the type of coefficients of powers of ω.

`ZOmega{T}` is defined as the alias
```julia
ZOmega{T} = CyclotomicRing{4, T} where {T <: Integer}
```
"""
const ZOmega{T} = CyclotomicRing{4, T} where {T <: Integer}

###
### Constructors
###

# This causes ambiguity with
# CyclotomicRing{M, CT1}(c::CyclotomicRing{M, CT2}) where {M, CT1 <: Number, CT2 <: Number}
function CyclotomicRing{4, T}(x::Number) where {T}
    z = zero(T)
    return CyclotomicRing{4, T}((T(x), z, z, z))
end

function CyclotomicRing{M, T}(x::Number) where {M, T}
    z = zero(T)
    val = T(x)
    return CyclotomicRing{M, T}(ntuple(i -> i == 1 ? val : z, Val(M)))
end

function CyclotomicRing()
    throw(MethodError(CyclotomicRing, ()))
end

function CyclotomicRing(coeff::T, coeffs::T...) where {T <: Number}
    CyclotomicRing((coeff, coeffs...))
end

function CyclotomicRing{M}(c1::Number, c2, coeffs...) where M
    cs = promote(c1, c2, coeffs...)
    CyclotomicRing{M, typeof(cs[1])}(cs)
end

function CyclotomicRing(c1::Dyadic, coeffs::Dyadic...)
    cs = promote(c1, coeffs...)
    CyclotomicRing{length(cs), typeof(cs[1])}(cs)
end

function CyclotomicRing{4}(a::Number,b,c,d)
    cs = promote(a,b,c,d)
    CyclotomicRing{4, typeof(cs[1])}(cs)
end

function DOmega(cs::NTuple{4, <:Dyadic})
    CyclotomicRing{4, typeof(first(cs))}(cs)
end

function DOmega(cs::NTuple{4, <:Integer})
    cs = map(Dyadic, cs)
    CyclotomicRing{4, typeof(first(cs))}(cs)
end

function DOmega(c1, c2, c3, c4)
    cs0 = (c1, c2, c3, c4)
    cs1 = promote(cs0...)
    cs = map(Dyadic, cs1)
    CyclotomicRing{4, Dyadic{typeof(cs[1].a), Int}}(cs)
end

function _mkomega(::Type{T}, a, b, c, d) where {T}
    coeffs = promote(a, b, c, d)
    coeffs1 = map(T, coeffs)
    V = typeof(coeffs1[1])
    CyclotomicRing{4, V}(coeffs)
end

DOmega{T}(a, b, c, d) where {T <: Integer} = _mkomega(Dyadic{T,Int}, a,b,c,d)

ZOmega(a::T, b::T, c::T, d::T) where {T <: Integer} = ZOmega{Int}(a, b, c, d)
ZOmega{T}(a, b, c, d) where {T <: Integer} = _mkomega(T, a,b,c,d)

ZOmega(a::Number) = ZOmega{Int}(a)

function ZOmega(c::CyclotomicRing{4, CT}) where {CT}
    coeffs = map(x -> convert(Integer, x), c.coeffs)
    ZOmega{typeof(first(coeffs))}(coeffs)
end

function DOmega(c::CyclotomicRing{4, CT}) where {CT}
    coeffs = map(x -> Dyadic(x), c.coeffs)
    T = typeof(first(coeffs).a)
    DOmega{T}(coeffs)
end

function DOmega(r::Rational)
    d = Dyadic(r)
    z = zero(d)
    T = typeof(d.a)
    DOmega{T}(d, z, z, z)
end

DOmega(a::T) where {T<:Integer} = DOmega{T}(a)

# Resolves ambiguity
function CyclotomicRing{4, CT1}(c::CyclotomicRing{4, CT2}) where {CT1 <: Number, CT2 <: Number}
    coeffs = convert.(CT1, c.coeffs)
    CyclotomicRing{4, CT1}(coeffs)
end

function CyclotomicRing{M, CT1}(c::CyclotomicRing{M, CT2}) where {M, CT1 <: Number, CT2 <: Number}
    coeffs = convert.(CT1, c.coeffs)
    CyclotomicRing{M, CT1}(coeffs)
end

###
### End constructors
###

###
### Show
###

function Base.show(io::IO, ::PRETTY, cr::CyclotomicRing)
    c = cr.coeffs
    n = length(c)
    showcount = 0
    for i in 1:n
        iszero_strong(c[i]) && continue
        showcount += 1
        if showcount > 1
            print(io, " + ")
        end
        if isone_strong(-c[i])
            print(io, "-")
        elseif !isone_strong(c[i])
            show(io, PRETTY(), c[i])
        end
        if isone_strong(i-1)
            print(io, "ω")
        else
            print(io, "ω", superscript(i-1))
        end
    end
    if showcount == 0
        show(io, PRETTY(), zero(first(c)))
        # TODO, use zero of first(c) somehow.
        # This will be more robust
#        show(io, MIME"text/plain", zero(first(c)))
    end
end

# This is just a desperate attempt to avoid unbound args error from Aqua
# function Base.show(io::IO, ::PRETTY, tup::NTuple{0, <:CyclotomicRing})
#     print(io, "()")
# end

#function Base.show(io::IO, ::PRETTY, tup::NTuple{<:Any, T}) where {T<:CyclotomicRing}

# function Base.show(io::IO, ::PRETTY, rings::T...) where {T<:CyclotomicRing}
#     print(io, "(")
#     for i in 1:length(rings)
#         if i > 1
#             print(io, ", ")
#         end
#         show(io, PRETTY(), rings[i])
#     end
#     print(io, ")")
# end

# Seems I have to do this piracy
# function Base.show(io::IO, m::MIME{Symbol("text/plain")}, t::Tuple{})
#     invoke(show, Tuple{IO,MIME{Symbol("text/plain")},Any}, stdout, m, t)
# end

#function Base.show(io::IO, ::PRETTY, rings::NTuple{<:Any, T}) where {T<:CyclotomicRing}
function Base.show(io::IO, ::PRETTY, rings::Tuple{T, Vararg{T, nMinus1}}) where {T<:CyclotomicRing, nMinus1}
    print(io, "(")
    for i in 1:length(rings)
        if i > 1
            print(io, ", ")
        end
        show(io, PRETTY(), rings[i])
    end
    print(io, ")")
end

###
### Accessors
###

"""
    coeffs(cyc::CyclotomicRing)

Return a `Tuple` of the coeffients of `cyc`

# Examples
```jldoctest
julia> coeffs(ZRoot2(1,2))
(1, 2)

julia> coeffs(DRoot2(3,4))
(3, 4)

julia> coeffs(DRoot2(1,Dyadic(3,2)))
(1, 3/2²)
```
"""
@inline coeffs(cyc::CyclotomicRing) = cyc.coeffs

"""
    Base.getindex(cyc::CyclotomicRing, n::Integer)

Return the coefficient of the `n`th-order term in the polynomial `cyc`.

The first term is the highest order term. The last term is of order zero.
So, `cyc[0]` returns the term of order zero. We are assuming that the primitive
root is an `N`th root of unity, with `N` even. A polynomial in such a primitive
root has only `N ÷ 2` coefficients.
"""
Base.getindex(cyc::CyclotomicRing, n::Integer) = cyc.coeffs[n+1]

Base.iterate(cyc::CyclotomicRing, i::Integer=1) = iterate(cyc.coeffs, i)
Base.eltype(cyc::CyclotomicRing{<:Any, T}) where {T} = T
Base.length(cyc::CyclotomicRing) = length(cyc.coeffs)
Base.lastindex(::CyclotomicRing{M}) where M = M - 1
Base.transpose(cyc::CyclotomicRing) = cyc


function Base.conj(cyc::DOmega{T}) where T
    (a, b, c, d) = cyc.coeffs
    # I think promotion not needed.
    newcoeff = promote(a, -d, -c, -b)
    DOmega{T}(newcoeff)
end

"""
    conj(cyc::CyclotomicRing{4})

Return the complex conjugate of `cyc`.
"""
function Base.conj(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    newcoeff = (a, -d, -c, -b)
    typeof(cyc)(newcoeff)
end

function ZRoot2(z::ZOmega)
    (a, b, c, d) = coeffs(z)
    iszero(c) && d == -b || throw(ArgumentError(lazy"Inexact error converting $z to ZRoot2."))
    ZRoot2(a, b)
end

function ZOmega(q::ZRoot2)
    ZOmega(q.a, q.b, zero(q.b), -q.b)
end

function promote_rule(::Type{<:ZRoot2}, ::Type{T}) where {T <: ZOmega}
    T
end

# This looks wrong
# In addition nothing seems to dispatch to this. Although I expect it should
function promote_rule(::Type{CyclotomicRing{<:Any, T}}, ::Type{V}) where {T, V <: Base.BitInteger}
    promote_type(promote_type(T, AbstractFloat), V)
end

# This works. Why we need the hard coded order, 4?
function promote_rule(::Type{CyclotomicRing{4, T}}, ::Type{V}) where {T, V <: Base.BitInteger}
    CyclotomicRing{4, promote_type(T, V)}
end

function promote_rule(::Type{T}, ::Type{<:RootOne{8}}) where {T}
    T
end
promote_rule(::Type{<:RootOne{8}}, ::Type{T}) where {T} = T


# Base.:(==)(a::ZRoot2, b::ZOmega) = ZOmega(a) == b
# Base.:(==)(a::ZOmega, b::ZRoot2) = b == a

Base.one(::Type{DOmega}) = one(DOmega{Int})
Base.one(::Type{ZOmega}) = one(ZOmega{Int})
Base.zero(::Type{DOmega}) = zero(DOmega{Int})
Base.zero(::Type{ZOmega}) = zero(ZOmega{Int})

Base.adjoint(cyc::CyclotomicRing) = conj(cyc)

LinearAlgebra.norm(cyc::CyclotomicRing) = abs(cyc)
Base.abs(cyc::CyclotomicRing) = sqrt(abs2(cyc))

function Base.abs2(cyc::CyclotomicRing)
    canonical(abs2(complex(cyc)))
end

    # def norm(self):
    #     return (self._a**2 + self._b**2 + self._c**2 + self._d**2) ** 2 - 2 * (
    #         self._a * self._b
    #         + self._b * self._c
    #         + self._c * self._d
    #         - self._d * self._a
    #     ) ** 2

function rnorm(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    x = a^2 + b^2 + c^2 + d^2
    y = a * b + b * c + c * d - a * d
    x*x - 2 * y * y
end

"""
    conj_root_two(cyc::CyclotomicRing{4})
    conj_root_two(cyc::DOmega)
    conj_root_two(cyc::ZOmega)

Maps `a + b√2` to `a - b√2` in all coefficients in `cyc`.

This implements √2-conjugation.
"""
function conj_root_two(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    newcoeff = (a, -b, c, -d)
    typeof(cyc)(newcoeff)
end

"""
    imaginary(::Type{DOmega{T}}) where {T}

The value of type `DOmega{T}` that represents the imaginary unit.
"""
function imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
    CyclotomicRing(z, z, o, z)
end

"""
    sqrt_imaginary(::Type{DOmega{T}}) where {T}

The value of type `DOmega{T}` that represents the principal square root of the imaginary unit.
"""
function sqrt_imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
    CyclotomicRing(z, o, z, z)
end

(t::Type{CyclotomicRing{4, T}})(::Type{RootImagT}) where {T} = sqrt_imaginary(t)

"""
    one_over_root_two(::Type{DOmega{T}}) where {T}

Return a value of `DOmega{T}` representing the reciprocal of the square root of two.

# Example
```jldoctest
julia> one_over_root_two(DOmega{Int})
1/2 ω + -1/2 ω³
```
"""
function one_over_root_two(::Type{CyclotomicRing{4, Dyadic{T, Int}}}) where {T}
    DT = Dyadic{T, Int}
    z = zero(DT)
    half = Dyadic(T(1), 1)
#    CyclotomicRing(-half, z, half, z)
    CyclotomicRing(z, half, z, -half)
end

"""
    root_two(::Type{DOmega{T}}) where {T}

Return a value of `DOmega{T}` representing the square root of two.

# Example
```jldoctest
julia> root_two(DOmega{Int})
 ω + - ω³
```
"""
function root_two(::Type{CyclotomicRing{4, Dyadic{T, Int}}}) where {T}
    DT = Dyadic{T, Int}
    z = zero(DT)
    o = one(DT)
#    CyclotomicRing(-half, z, half, z)
    CyclotomicRing(z, o, z, -o)
end

function Base.real(cyc::ZOmega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    DRoot2(Dyadic(a, 0), Dyadic(b-d, 1))
end

function Base.real(cyc::DOmega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    DRoot2(a, mul_half(b - d))
end

function Base.imag(cyc::ZOmega{T}) where T
    (a, b , c, d) = coeffs(cyc)
    DRoot2(Dyadic(c, 0), Dyadic(b+d, 1))
end

function Base.imag(cyc::DOmega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    DRoot2(c, mul_half(b + d))
end

#function Base.isreal(cyc::DOmega{<:Integer})
function Base.isreal(cyc::CyclotomicRing{4}, approx::AbstractApprox=Equal())
    (a, b , c, d) = coeffs(cyc)
    iszero(c) && iszero(b + d)
end

DOmega(zz::Complex{<:DRoot2}) = CyclotomicRing{4}(zz)
CyclotomicRing(zz::Complex{<:DRoot2}) = CyclotomicRing{4}(zz)
function CyclotomicRing{4}(zz::Complex{<:DRoot2})
    (r, i) = (real(zz), imag(zz))
    (w, x) = (r.a, r.b)
    (y, z) = (i.a, i.b)
    a = w
    c = y
    b = canonical(x + z)
    d = canonical(z - x)
    return DOmega((a,b,c,d))
end

DOmega(zz::DRoot2) = CyclotomicRing{4}(zz)
CyclotomicRing(zz::DRoot2) = CyclotomicRing{4}(zz)
function CyclotomicRing{4}(zz::DRoot2)
    (w, x) = (zz.a, zz.b)
    a = w
    c = zero(typeof(w))
    b = canonical(x)
    d = canonical(- x)
    return DOmega((a,b,c,d))
end

function RootOne{8}(cyc::CyclotomicRing{4})
    cs = coeffs(cyc)
    if (count(x -> isone(abs2(x)), cs) == 1) &&
        (count(iszero, cs) == 3)
        ind = findfirst(!=(0), cs)
        sgn = Int(cs[ind])
        rt = RootOne{8}(ind  - 1)
        return sgn == 1 ? rt : -rt
    end
end

function canonical(c::CyclotomicRing)
    CyclotomicRing(map(canonical,  c.coeffs))
end

function Base.:+(c1::CyclotomicRing{M, CoeffT}, c2::CyclotomicRing{M, CoeffT}) where {M, CoeffT}
    CyclotomicRing{M, CoeffT}(c1.coeffs .+ c2.coeffs)
end

function Base.:-(c1::CyclotomicRing{M}, c2::CyclotomicRing{M}) where {M}
    coeffs = promote((c1.coeffs .- c2.coeffs)...)
    T = typeof(first(coeffs))
    CyclotomicRing{M, T}(coeffs)
end

function Base.:+(c1::CyclotomicRing{M}, c2::CyclotomicRing{M}) where {M}
    coeffs = promote((c1.coeffs .+ c2.coeffs)...)
    T = typeof(first(coeffs))
    CyclotomicRing{M, T}(coeffs)
end

function Base.:-(c::CyclotomicRing)
    CyclotomicRing(.- c.coeffs)
end

function Base.:*(c1::CyclotomicRing{4}, c2::CyclotomicRing{4})
    (a1, b1, c1, d1) = c1.coeffs
    (a2, b2, c2, d2) = c2.coeffs
    coeffs = (
        a1*a2 - c1*c2 - d1*b2 - b1*d2,

        a1*b2 + b1*a2 - c1*d2 - d1*c2,

        a1*c2 + b1*b2 + c1*a2 - d1*d2,

        a1*d2 + b1*c2 + c1*b2 + d1*a2,
    )
    CyclotomicRing{4, typeof(first(coeffs))}(coeffs)
end

function Base.:*(r::Rational, c::CyclotomicRing)
    Dyadic(r) * c
end

function Base.:*(r::Real, c::CyclotomicRing)
    # These two are same speed.
    #    cs = map(x -> r * x, c.coeffs)
    cs = c.coeffs .* r
    # Use NTuple to avoid big perf hit
    # cs = Tuple(r * x for x in c.coeffs) slow
    CyclotomicRing(cs)
end

@inline Base.:*(c::CyclotomicRing, r::Real) = r * c

# Wrote this all out explicitly in case there was a broadcasting inefficiency or s.t
# But no!
@inline function Base.:*(c::CyclotomicRing{4}, z::Complex)
    (rz, iz) = (real(z), imag(z))
    (a, b, c, d) = c.coeffs
    # Contribution from real part
    (a2, b2, c2, d2) = (rz * a, rz * b, rz * c, rz * d)
    # Contribution from imag part
    (a1, b1, c1, d1) = (iz * a, iz * b, iz * c, iz * d)

    # i permutes and also changes some signs. Not quite automorphism
#    (ai, bi, ci, di) = (c1, d1, -a1, -b1)
    (ai, bi, ci, di) = (-c1, -d1, a1, b1)

    # Add the two contributions
    cs = (a2 + ai, b2 + bi, c2 + ci, d2 + di)
    CyclotomicRing(cs)
end

@inline Base.:*(z::Complex, c::CyclotomicRing) = c * z

function Base.:*(x::DRoot2, c::CyclotomicRing)
    error("not implemented")
end

function Base.:+(r::RootOne{8}, s::RootOne{8})
    DOmega(r) + DOmega(s)
end

function Base.:-(r::RootOne{8}, s::RootOne{8})
    ZOmega(r) - ZOmega(s)
end

function Base.:*(x::Integer, r::RootOne{8})
    x * ZOmega(r)
end

function Base.:*(x::Dyadic, r::RootOne{8})
    x * DOmega(r)
end

function Base.:*(x::Rational, r::RootOne{8})
    Dyadic(x) * DOmega(r)
end

function Base.:*(r::RootOne{8}, cyc::CyclotomicRing{4})
    k = r.k
    (a, b, c, d) = cyc.coeffs
    coeffs =
        if k == 0
            (a,b,c,d)     # 1
        elseif k == 1
            (-d,a,b,c)    # √𝕚
        elseif k == 2
            (-c,-d,a,b)   # 𝕚
        elseif k == 3
            (-b,-c,-d,a)  # 𝕚√𝕚
        elseif k == 4
            (-a,-b,-c,-d) # -1
        elseif k == 5
            (d,-a,-b,-c)  # -√𝕚
        elseif k == 6
            (c,d,-a,-b)   # -𝕚
        elseif k == 7
            (b,c,d,-a)    # -𝕚√𝕚
        end
    return typeof(cyc)(coeffs)
end

Base.:*(c::CyclotomicRing{4}, r::RootOne{8}) = r * c

Base.:/(c::CyclotomicRing{4}, r::RootOne{8}) = c * inv(r)

function Base.:*(::ImagT, cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    typeof(cyc)(-c,-d,a,b)
end

function Base.:*(::RootImagT, cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    typeof(cyc)(-d,a,b,c)
end

# This expression is invariant wrt reversing
# the order of coefficients.
"""
    mul_root_two(cyc::CyclotomicRing{4}, n::Integer=1)

Multiply `cyc` by `√2ⁿ`
"""
function mul_root_two(cyc::CyclotomicRing{4})
    (a,b,c,d) = cyc.coeffs
#    print("0")
    coeffs = (b-d, c+a, b+d, c-a)
    typeof(cyc)(coeffs)
end

function Base.:*(pow::Pow{RootTwoT}, cyc::CyclotomicRing{4})
    mul_root_two(cyc, pow.n)
end

function mul_root_two(cyc::CyclotomicRing{4}, n::Integer)
    n == 0 && return cyc
    if n > 0
        # We might want something more efficient, less overflow prone, etc.
        cyc1 = mul_two(cyc, div(n,2))
        iseven(n) && return cyc1
        return mul_root_two(cyc1)
    else
        n = -n
        iseven(n) && return mul_half(cyc, div(n, 2))
        return mul_root_two(mul_half(cyc, div(n, 2) + 1))
    end
end

# This converts to float. Not what we want, I think.
# We need to not change the type of the input
"""
    mul_half(cyc::CyclotomicRing{4}, n::Integer=1)

Multiply `cyc` by `n` factors of the reciprocal of two.

`n` may be positive, negative, or zero.
"""
function mul_half(cyc::CyclotomicRing{4, <:Dyadic}, n::Integer=1)
    n == 0 && return cyc
    CyclotomicRing(map(x -> mul_half(x, n), cyc.coeffs))
end

# TODO: This should check that each coeff is divisible by 2^n
# Intentionally convert ZOmega to DOmega
function mul_half(cyc::CyclotomicRing{4, <: Integer}, n::Integer=1)
    CyclotomicRing(map(x -> Dyadic(x, n), coeffs(cyc)))
    # newc = map(x -> div(x, 2^n), coeffs(cyc))
    # @show newc
    # CyclotomicRing(newc)
end

# BAD BAD
# function mul_half(cyc::CyclotomicRing{4, <: Integer})
#     CyclotomicRing(map(x -> div(x, 2), coeffs(cyc)))
# end

function div_coefficients(cyc::CyclotomicRing{4, <: Integer}, fac)
    typeof(cyc)(map(x -> div(x, fac), coeffs(cyc)))
end

# BAD!
function div_half(cyc::CyclotomicRing{4, <: Integer})
    CyclotomicRing(map(x -> div(x, 2), coeffs(cyc)))
end

@inline function Base.:*(::InvTwoT, cyc::CyclotomicRing{4}, n::Integer=1)
    mul_half(cyc)
end

function Base.:*(p::Pow{InvTwoT}, cyc::CyclotomicRing{4, <: Integer}, n::Integer=1)
    mul_half(cyc, p.n)
end

@inline function Base.:*(::TwoT, cyc::CyclotomicRing{4})
    mul_two(cyc)
end

Base.:*(cyc::CyclotomicRing{4}, ::TwoT) = Two * cyc

@inline function Base.:*(powtwo::Pow{TwoT}, cyc::CyclotomicRing{4})
    mul_two(cyc, powtwo.n)
end

@inline Base.:*(cyc::CyclotomicRing{4}, powtwo::Pow{TwoT}) = powtwo * cyc

@inline function mul_two(cyc::CyclotomicRing{4}, n::Integer=1)
    n == 0 && return cyc
    CyclotomicRing(map(x -> mul_two(x, n), cyc.coeffs))
end

# function Base.:*(pow::Pow, cyc::CyclotomicRing{4})
#     CyclotomicRing(map(x -> pow * x, cyc.coeffs))
# end

# sqrt(2) == omega - omega^3
function Base.:*(::RootTwoT, cyc::CyclotomicRing{4})
    (a,b,c,d) = cyc.coeffs
    coeffs = (b-d, c+a, b+d, c-a)
    CyclotomicRing(coeffs)
end

Base.:*(::InvRootTwoT, cyc::CyclotomicRing{4}) = mul_one_over_root_two(cyc)

function Base.:*(p::Pow{InvRootTwoT}, cyc::CyclotomicRing{4})
    mul_one_over_root_two(cyc, p.n)
end

function mul_one_over_root_two(cyc::CyclotomicRing{4}, n::Integer)
    if iseven(n)
        mul_half(cyc, div(n, 2))
    else
        mul_one_over_root_two(mul_half(cyc, div(n, 2)))
    end
end

"""
    mul_one_over_root_two(cyc::CyclotomicRing{4})

Multiply `cyc` by the reciprocal of the square root of two.
"""
function mul_one_over_root_two(cyc::CyclotomicRing{4})
    mul_root_two(mul_half(cyc))
end

@inline Base.:(==)(c1::CyclotomicRing{N}, c2::CyclotomicRing{N}) where {N} = c1.coeffs == c2.coeffs

function Base.:(==)(c1::CyclotomicRing, c2::CyclotomicRing)
    throw(ArgumentError(lazy"Unsupported comparison"))
end

@inline function Base.:(==)(cyc::CyclotomicRing{N}, r::RootOne{M}) where {N, M}
    2*N == M || return throw(ArgumentError(lazy"Unsupported arguments"))
    r.k > N && return -cyc == -r
    isone(sum(isone, coeffs(cyc))) || return false # Exactly one coeff is one
    n = findfirst(isone, coeffs(cyc)) # Which coeff is one?
    n == r.k + 1
end

function Base.:^(c::CyclotomicRing, n::Integer)
    n == 0 && return one(c)
    n == 1 && return c
    n == 2 && return c * c
    return Base.power_by_squaring(c, n)
end

function Base.one(::Type{CyclotomicRing{D, CoeffT}}) where {D, CoeffT}
    CyclotomicRing{D, CoeffT}(one(CoeffT), ntuple(x -> zero(CoeffT), D-1)...,)
end

Base.isone(cyc::CyclotomicRing) = canonical(cyc) == one(cyc)

Base.one(::Type{CyclotomicRing{D}}) where {D} = one(CyclotomicRing{D, Int})

Base.one(c::CyclotomicRing) = one(typeof(c))
Base.zero(::Type{CyclotomicRing{D, CoeffT}}) where {D, CoeffT} = CyclotomicRing(ntuple(x -> zero(CoeffT), D)...,)
Base.zero(::Type{CyclotomicRing{D}}) where {D} = zero(CyclotomicRing{D, Int})
Base.zero(cr::CyclotomicRing) = zero(typeof(cr))

function Base.convert(::Type{CyclotomicRing{M, CT1}}, c::CyclotomicRing{M, CT2}) where {M, CT1, CT2}
    CyclotomicRing{M, CT1}(c)
end

# This is a stopgap fix for a bug.
# The more generic methods are introducing a Float64 somewhere in the
# conversion. They should be fixed. But this method explicitly catches the case of BigFloat
Base.float(cyc::CyclotomicRing{<:Any, BigFloat}) = big(cyc)

Base.float(cyc::CyclotomicRing{4, BigFloat}) = Complex{BigFloat}(cyc)

function Base.big(c::CyclotomicRing)
    Complex{BigFloat}(c)
end

function Base.Complex(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    r8 = RootOne{8}
    float(a) + float(b) * float(r8(1)) + float(c) * float(r8(2)) +
        float(d) * float(r8(3))
end

function Base.Complex(cyc::CyclotomicRing{4,T}) where {T <: Union{Integer, Dyadic}}
    Complex(real(cyc), imag(cyc))
end

Base.float(cyc::CyclotomicRing{4}) = Complex(float(real(cyc)), float(imag(cyc)))
Base.complex(cyc::CyclotomicRing) = Complex(cyc)

function Base.Complex(cyc::CyclotomicRing{D, T}) where {D, T}
    Tc = float(T)
    _convert_cyc(cyc, Val(D), Complex{Tc})
end

function Base.Complex{Tc}(cyc::CyclotomicRing{D}) where {D, Tc}
    _convert_cyc(cyc, Val(D), Complex{Tc})
end

function (::Type{T})(cyc::CyclotomicRing{D}) where {T<:Real, D}
    (f, rest...) = coeffs(cyc)
    all(iszero, rest) || throw(ArgumentError(lazy"Inexact error converting $cyc to $T."))
    T(f)
end

function Complex{T}(cyc::CyclotomicRing{D}) where {T<:Integer, D}
    (f, rest...) = coeffs(cyc)
    all(iszero, rest) || throw(ArgumentError(lazy"Inexact error converting $cyc to Complex{$T}."))
    Complex{T}(f)
end

function Complex{T}(cyc::CyclotomicRing{4}) where {T<:Integer}
    (a, b, c, d) = coeffs(cyc)
    if all(iszero, (b, c, d))
        return Complex{T}(a)
    end
    if all(iszero, (a, b, d))
        return Complex{T}(zero(T), T(c))
    end
    if all(iszero, (b, d))
        return Complex{T}(a) + Complex{T}(zero(T), T(c))
    end
    throw(ArgumentError(lazy"Inexact error converting $cyc to Complex{$T}."))
end

Base.float(cyc::CyclotomicRing) = complex(cyc)

# Faster than using any kind of iterative or reduce scheme
# This is only slightly slower 2.7 vs 3.5 ns than the explicit
# expression above for D = 4
@inline function _convert_cyc(c, ::Val{D}, func) where {D}
    _add_term(func, 0, Val(D), coeffs(c)...)
end

@inline function _add_term(func, i, ::Val{D}, c) where {D}
    func(RootOne{2*D}(i)) * func(c)
end

@inline function _add_term(func, i, ::Val{D}, c, _coeffs...) where {D}
    func(RootOne{2*D}(i)) * func(c) + _add_term(func, i+1, Val(D), _coeffs...)
end

# Slightly faster than the general method: 3.6ns vs 2.6ns
function Base.Complex{Tc}(c::CyclotomicRing{4}) where Tc
    (a, b, c, d) = c.coeffs
    r8 = RootOne{8}
    T = Complex{Tc}
    T(a) + T(b) * T(r8(1)) + T(c) * T(r8(2)) +
        T(d) * T(r8(3))
end

# Huh. I did this twice

ZOmega(r::RootOne{8}) = CyclotomicRing{4, Int}(r)
DOmega(r::RootOne{8}) = DOmega{Int}(r)

function CyclotomicRing{4, T}(r::RootOne{8}) where {T}
    val = r.k == 4 ? -1 : sign(4 - r.k)
    pos = mod(r.k, 4)
    coeffs =
        if pos == 0
            (val,0,0,0)
        elseif pos == 1
            (0,val,0,0)
        elseif pos == 2
            (0,0,val,0)
        else
            (0,0,0,val)
        end
    #    println(T)
    cs = map(T, coeffs)
    CyclotomicRing{4, typeof(cs[1])}(cs)
end

# function CyclotomicRing{4}(r::RootOne{8})
#     k = r.k
#     coeffs = if k == 0
#         (1, 0, 0, 0)
#     elseif k == 1
#         (0, 1, 0, 0)
#     elseif k == 2
#         (0, 0, 1, 0)
#     elseif k == 3
#         (0, 0, 0, 1)
#     elseif k == 4
#         (-1, 0, 0, 0)
#     elseif k == 5
#         (0, -1, 0, 0)
#     elseif k == 6
#         (0, 0, -1, 0)
#     elseif k == 7
#         (0, 0, 0, -1)
#     end
#     CyclotomicRing{4}(coeffs...)
# end

# With respect to base sqrt(2)
# Eg. lde( (omega - omega^3)/2) is one, not two
function least_denominator_exponent(cyc::CyclotomicRing{4, <:Dyadic})
    ccyc = canonical(cyc)
    ks = map(x -> x.k, ccyc.coeffs)
    (a, b, c, d) = map(x -> x.a, ccyc.coeffs)
    k = maximum(ks)
    (ka, kb, kc, kd) = ks
    z = zero(a)
    ap = k == ka ? a : z
    bp = k == kb ? b : z
    cp = k == kc ? c : z
    dp = k == kd ? d : z
    if k > 0 && iseven(ap - cp) && iseven(bp - dp)
        2 * k - 1
    else
        2 * k
    end
end

"""
    struct DOmegaA{T, KT} <: Number

An element of `𝔻[ω] = ℤ[1/√2, i]` represented as an element of `ℤ[ω]` together with
 a factor of a power of `1/√2`.

# Examples
```jldoctest
julia> z = DOmega(3//4,1,1,1//2)
3/2²ω⁰ + ω + ω² + 1/2ω³

julia> za = DOmegaA(z)
(3ω⁰ + 4ω + 4ω² + 2ω³) / √2⁴

julia> DOmega(za)
3/2²ω⁰ + ω + ω² + 1/2ω³

julia> DOmega(za) === z
true
```
"""
struct DOmegaA{T, KT} <: Number
    # No, this makes no sense
    # function DOmegaA(z::ZOmega{T}, k::KT) where {T, KT}
    #     iszero(least_denominator_exponent(z)) || throw(ArgumentError(lazy"Bad z constructing DOmegaA "))
    #     new{T, KT}(z, k)
    # end

    z::ZOmega{T}
    k::KT
end

DOmegaA(z::ZOmega) = DOmegaA(z, 0)

# This should be internal, so that you can't consruct an
# Illegal one.
@inline function DOmegaA(z::DOmega)
    lde = least_denominator_exponent(z)
    DOmegaA(ZOmega(mul_root_two(z, lde)), lde)
end

function least_denominator_exponent(z::DOmegaA)
    return z.k
end

function Base.show(io::IO, ::PRETTY, x::DOmegaA)
    print(io, "(")
    show(io, PRETTY(), x.z)
    print(io, ")")
    print(io, " / √2", superscript(x.k))
end

function DOmega(z::DOmegaA)
    DOmega(map(x -> canonical(Dyadic(x, div(z.k, 2))), coeffs(z.z)))
end

canonical(z::DOmegaA) = DOmegaA(canonical(z.z), z.k)

function Base.:(==)(c1::DOmegaA, c2::DOmegaA)
    (d1, d2) = (canonical(c1), canonical(c2))
    d1.z == d2.z && d1.k == d2.k
end

Base.:(==)(c1::DOmegaA, c2::DOmega) = c1 == DOmegaA(c2)

function ZOmega(z::DOmegaA)
    z = canonical(z)
    iszero(z.k) || throw(ArgumentError(lazy"Inexact error converting $z to ZOmega."))
    z.z
end

Base.conj(z::DOmegaA) = DOmegaA(conj(z.z), z.k)
Base.one(z::DOmegaA) = DOmegaA(one(z.z), zero(z.k))
Base.zero(z::DOmegaA) = DOmegaA(zero(z.z), zero(z.k))
Base.transpose(z::DOmegaA) = z

Base.real(z::DOmegaA) = mul_half(real(z.z), z.k >> 1)
Base.imag(z::DOmegaA) = mul_half(imag(z.z), z.k >> 1)
Base.Complex(z::DOmegaA) = Complex(real(z), imag(z))
Base.complex(z::DOmegaA) = Complex(z)

DOmegaA(zz::Complex{<:DRoot2}) = DOmegaA(DOmega(zz))

end # module CyclotomicRings
