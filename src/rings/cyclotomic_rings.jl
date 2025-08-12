module CyclotomicRings

import LinearAlgebra
import Base: convert, zero, one, promote_rule
import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two, root_two, coeffs,
    mul_root_two, mul_one_over_root_two, mul_half, conj_root_two, mul_two
import ..Utils: superscript, iszero_strong, isone_strong, PRETTY
import ..RootOnes: RootOne

import ..Dyadics: Dyadic
import ..QuadraticRings: Droot2
import ..Singletons: InvTwo, InvTwoT,
    RootTwo, RootTwoT,
    InvRootTwo, InvRootTwoT,
    Imag, ImagT,
    RootImag, RootImagT,
    Pow

export CyclotomicRing, Zomega, Domega

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
```
const Zomega{T} = CyclotomicRing{4, T} where {T <: Integer}
```

The ring `𝔻[ω] = ℤ[1/√2, i]`:
```
const Domega{T} = CyclotomicRing{4, Dyadic{T, Int}}
```
"""
struct CyclotomicRing{M, CoeffT}
    coeffs::NTuple{M, CoeffT}
end

"""
    coeffs(cyc::CyclotomicRing)

Return a `Tuple` of the coeffients of `cyc`

# Examples
```jldoctest
julia> coeffs(Zroot2(1,2))
(1, 2)

julia> coeffs(Droot2(3,4))
(3, 4)

julia> coeffs(Droot2(1,Dyadic(3,2)))
(1, 3/2²)
```
"""
@inline coeffs(cyc::CyclotomicRing) = cyc.coeffs

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
            print(io, " ω")
        else
            print(io, " ω", superscript(i-1))
        end
    end
    if showcount == 0
        show(io, PRETTY(), zero(first(c)))
        # TODO, use zero of first(c) somehow.
        # This will be more robust
#        show(io, MIME"text/plain", zero(first(c)))
    end
end

function Base.show(io::IO, ::PRETTY, tup::NTuple{<:Any, T}) where {T<:CyclotomicRing}
    print(io, "(")
    for i in 1:length(tup)
        if i > 1
            print(io, ", ")
        end
        show(io, PRETTY(), tup[i])
    end
    print(io, ")")
end


# FIX THIS! We want to agree with the rest of the world.
"""
    Base.getindex(cyc::CyclotomicRing, n::Integer)

Return the coefficient of the `n`th-order term in the polynomial `cyc`.

The first term is the highest order term. The last term is of order zero.
So, `cyc[0]` returns the term of order zero. We are assuming that the primitive
root is an `N`th root of unity, with `N` even. A polynomial in such a primitive
root has only `N ÷ 2` coefficients.
"""
function Base.getindex(cyc::CyclotomicRing, n::Integer)
    cyc.coeffs[n+1]
end

Base.iterate(cyc::CyclotomicRing, i::Integer=1) = iterate(cyc.coeffs, i)
Base.eltype(cyc::CyclotomicRing{<:Any, T}) where {T} = T
Base.length(cyc::CyclotomicRing) = length(cyc.coeffs)
Base.lastindex(::CyclotomicRing{M}) where M = M - 1
Base.transpose(cyc::CyclotomicRing) = cyc

# This looks wrong
function promote_rule(::Type{CyclotomicRing{<:Any, T}}, ::Type{V}) where {T, V <: Base.BitInteger}
    promote_type(promote_type(T, AbstractFloat), V)
end

# The type of the exponent of (1/2) is hardcoded to Int.
# We probably only need one type for this field.
"""
    Domega{T}

Represents the ring `𝔻[ω] = ℤ[1/√2, i]`.

The type `T<:Integer` is the type of the numerator in the dyadic fractions.

Here, `𝔻 = ℤ[½]` is the ring of dyadic fractions, implemented by `Dyadic`.

`Domega{T}` is defined as the alias
```
Domega{T} = CyclotomicRing{4, Dyadic{T, Int}}
```
"""
const Domega{T} = CyclotomicRing{4, Dyadic{T, Int}}

"""
    Zomega{T <: Integer}

Represents the ring `ℤ[ω]`.

The type `T<:Integer` is the type of coefficients of powers of ω.

`Zomega{T}` is defined as the alias
```
Zomega{T} = CyclotomicRing{4, T} where {T <: Integer}
```
"""
const Zomega{T} = CyclotomicRing{4, T} where {T <: Integer}

function CyclotomicRing{4, T}(x::Number) where {T}
    z = zero(T)
    return CyclotomicRing{4, T}((T(x), z, z, z))
end

function Base.conj(cyc::Domega{T}) where T
    (a, b, c, d) = cyc.coeffs
    # I think promotion not needed.
    newcoeff = promote(a, -d, -c, -b)
    Domega{T}(newcoeff)
end

Base.one(::Type{Domega}) = one(Domega{Int})
Base.one(::Type{Zomega}) = one(Zomega{Int})
Base.zero(::Type{Domega}) = zero(Domega{Int})
Base.zero(::Type{Zomega}) = zero(Zomega{Int})

"""
    conj(cyc::CyclotomicRing{4})

Return the complex conjugate of `cyc`.
"""
function Base.conj(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    newcoeff = (a, -d, -c, -b)
    typeof(cyc)(newcoeff)
end

Base.adjoint(cyc::CyclotomicRing) = conj(cyc)

LinearAlgebra.norm(cyc::CyclotomicRing) = abs(cyc)
Base.abs(cyc::CyclotomicRing) = sqrt(abs2(cyc))

function Base.abs2(cyc::CyclotomicRing)
    (a, b, c, d) = cyc.coeffs
    a1 = a^2 + b^2 + c^2 + d^2
    b1 = a*(b-d) + c*(b+d)
    c1 = zero(a)
    d1 = -b1
    CyclotomicRing((a1,b1,c1,d1))
end

"""
    conj_root_two(cyc::CyclotomicRing{4})
    conj_root_two(cyc::Domega)
    conj_root_two(cyc::Zomega)

Maps `a + b√2` to `a - b√2` in all coefficients in `cyc`.

This implements √2-conjugation.
"""
function conj_root_two(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    newcoeff = (a, -b, c, -d)
    typeof(cyc)(newcoeff)
end

"""
    imaginary(::Type{Domega{T}}) where {T}

The value of type `Domega{T}` that represents the imaginary unit.
"""
function imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
#    CyclotomicRing(z, o, z, z)
    CyclotomicRing(z, z, o, z)
end


"""
    sqrt_imaginary(::Type{Domega{T}}) where {T}

The value of type `Domega{T}` that represents the principal square root of the imaginary unit.
"""
function sqrt_imaginary(::Type{CyclotomicRing{4, T}}) where {T}
    z = zero(T)
    o = one(T)
#    CyclotomicRing(z, z, o, z)
    CyclotomicRing(z, o, z, z)
end

(t::Type{CyclotomicRing{4, T}})(::Type{RootImagT}) where {T} = sqrt_imaginary(t)

"""
    one_over_root_two(::Type{Domega{T}}) where {T}

Return a value of `Domega{T}` representing the reciprocal of the square root of two.

# Example
```jldoctest
julia> one_over_root_two(Domega{Int})
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
    root_two(::Type{Domega{T}}) where {T}

Return a value of `Domega{T}` representing the square root of two.

# Example
```jldoctest
julia> root_two(Domega{Int})
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

CyclotomicRing(coeffs::T...) where T = CyclotomicRing{length(coeffs), T}(coeffs)

function CyclotomicRing(c1, coeffs...)
    cs = promote(c1, coeffs...)
    CyclotomicRing(cs)
end

function CyclotomicRing{M}(c1, coeffs...) where M
    cs = promote(c1, coeffs...)
    CyclotomicRing{M, typeof(cs[1])}(cs)
end

function CyclotomicRing{4}(a,b,c,d)
    cs = promote(a,b,c,d)
    CyclotomicRing{4, typeof(cs[1])}(cs)
end

function Domega(c1, c2, c3, c4)
    c1a = Dyadic{typeof(c1), Int}(c1)
    cs = promote(c1a, c2, c3, c4)
    CyclotomicRing(cs)
end


# Not using this. So comment out.
# function CyclotomicRing{4, T}(a, b, c, d) where T
#     cs = (T(a),T(b),T(c),T(d))
#     CyclotomicRing(cs)
# end

function _mkomega(::Type{T}, a, b, c, d) where {T}
    coeffs = promote(a, b, c, d)
    coeffs1 = map(T, coeffs)
    V = typeof(coeffs1[1])
    CyclotomicRing{4, V}(coeffs)
end

#Domega(a, b, c, d) = Domega{Int}(a,b,c,d)

Domega{T}(a, b, c, d) where {T <: Integer} = _mkomega(Dyadic{T,Int}, a,b,c,d)

# Zomega(a, b, c, d) = Zomega{Int}(a,b,c,d)
Zomega(a::T, b::T, c::T, d::T) where {T <: Integer} = Zomega{Int}(a, b, c, d)
Zomega{T}(a, b, c, d) where {T <: Integer} = _mkomega(T, a,b,c,d)

Zomega(a::Number) = Zomega{Int}(a)

function Zomega(c::CyclotomicRing{4, CT}) where {CT}
    coeffs = map(x -> convert(Integer, x), c.coeffs)
    Zomega{typeof(first(coeffs))}(coeffs)
end

function Domega(c::CyclotomicRing{4, CT}) where {CT}
    coeffs = map(x -> Dyadic(x), c.coeffs)
    T = typeof(first(coeffs).a)
    Domega{T}(coeffs)
end

function Domega(r::Rational)
    d = Dyadic(r)
    z = zero(d)
    T = typeof(d.a)
    Domega{T}(d, z, z, z)
end

Domega(a::T) where {T<:Integer} = Domega{T}(a)
#Domega(a::Number) = Domega{Int}(a)

function Base.real(cyc::Zomega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    Droot2(Dyadic(a, 0), Dyadic(b-d, 1))
end

function Base.real(cyc::Domega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    Droot2(a, mul_half(b - d))
end

function Base.imag(cyc::Zomega{T}) where T
    (a, b , c, d) = coeffs(cyc)
    Droot2(Dyadic(c, 0), Dyadic(b+d, 1))
end

function Base.imag(cyc::Domega{<:Integer})
    (a, b , c, d) = coeffs(cyc)
    Droot2(c, mul_half(b + d))
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


# function Base.:-(c1::CyclotomicRing{M, CoeffT}, c2::CyclotomicRing{M, CoeffT}) where {M, CoeffT}
#     coeffs = promote((c1.coeffs .- c2.coeffs)...)
#     T = typeof(first(coeffs))
#     CyclotomicRing{M, T}(coeffs)
# end


function Base.:-(c::CyclotomicRing)
    CyclotomicRing(.- c.coeffs)
end

#function Base.:*(c1::CyclotomicRing{4, T1}, c2::CyclotomicRing{4, T2}) where {T1, T2}
function Base.:*(c1::CyclotomicRing{4}, c2::CyclotomicRing{4})
    (a1, b1, c1, d1) = c1.coeffs
    (a2, b2, c2, d2) = c2.coeffs
    coeffs = (
        a1*a2 - c1*c2 - d1*b2 - b1*d2,

        a1*b2 + b1*a2 - c1*d2 - d1*c2,

        a1*c2 + b1*b2 + c1*a2 - d1*d2,

        a1*d2 + b1*c2 + c1*b2 + d1*a2,
    )
    CyclotomicRing(coeffs)
end

function Base.:*(r::Real, c::CyclotomicRing)
    # These two are same speed.
    #    cs = map(x -> r * x, c.coeffs)
    cs = c.coeffs .* r

    # cs = Tuple(r * x for x in c.coeffs) slow
    CyclotomicRing(cs)
end

# Specialized method does not help
# function Base.:*(r::Real, c::CyclotomicRing{4})
#     (a, b, c, d) = c.coeffs
#     CyclotomicRing(r * a, r * b, r * c, r * d)
# end


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

function Base.:*(x::Droot2, c::CyclotomicRing)
    error("not implemented")
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
    mul_root_two(cyc::CyclotomicRing{4})

Multiply `cyc` by the square root of two.
"""
function mul_root_two(cyc::CyclotomicRing{4})
    (a,b,c,d) = cyc.coeffs
    coeffs = (b-d, c+a, b+d, c-a)
    CyclotomicRing(coeffs)
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

"""
    mul_half(cyc::CyclotomicRing{4}, n::Integer=1)

Multiply `cyc` by `n` factors of the reciprocal of two.

`n` may be positive, negative, or zero.
"""
function mul_half(cyc::CyclotomicRing{4}, n::Integer=1)
    n == 0 && return cyc
    CyclotomicRing(map(x -> mul_half(x, n), cyc.coeffs))
end

function mul_two(cyc::CyclotomicRing{4}, n::Integer=1)
    n == 0 && return cyc
    CyclotomicRing(map(x -> mul_two(x, n), cyc.coeffs))
end

function Base.:*(::InvTwoT, cyc::CyclotomicRing)
    CyclotomicRing(map(x -> InvTwo * x, cyc.coeffs))
end

function Base.:*(pow::Pow, cyc::CyclotomicRing{4})
    CyclotomicRing(map(x -> pow * x, cyc.coeffs))
end

function Base.:*(::RootTwoT, cyc::CyclotomicRing{4})
    (a,b,c,d) = cyc.coeffs
    coeffs = (b-d, c+a, b+d, c-a)
    CyclotomicRing(coeffs)
end

Base.:*(::InvRootTwoT, cyc::CyclotomicRing{4}) = RootTwo * (InvTwo * cyc)

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
    isone(sum(isone, cyc)) || return false # Exactly one coeff is one
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

# function Base.convert(::Type{T}, c::CyclotomicRing) where {T}
#     CyclotomicRing(map(x -> convert(T, x), c.coeffs))
# end

# We want less of convert
# function Base.convert(::Type{Complex{Tc}}, c::CyclotomicRing) where Tc
#     Complex{Tc}(c)
# end

function CyclotomicRing{M, CT1}(c::CyclotomicRing{M, CT2}) where {M, CT1, CT2}
    coeffs = convert.(CT1, c.coeffs)
    CyclotomicRing{M, CT1}(coeffs)
end

# This is the right idea. But and it produces the right Complex.
# But somehow errors are intrduced when converting to BigFloat.
# The problem could be fixed.
# function Base.Complex(cyc::CyclotomicRing{4})
#     Complex(real(cyc), imag(cyc))
# end

function Base.Complex(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
    r8 = RootOne{8}
    float(a) + float(b) * float(r8(1)) + float(c) * float(r8(2)) +
        float(d) * float(r8(3))
end

function Base.Complex(cyc::CyclotomicRing{4,BigInt})
    (a, b, c, d) = cyc.coeffs
    cfunc = big
    r8 = RootOne{8}
    float(a) + float(b) * big(r8(1)) + float(c) * big(r8(2)) +
        float(d) * big(r8(3))
end

#Base.Complex(cyc::CyclotomicRing{4}) = convert(Complex, cyc)

# This gives incorrect result for BigInt
# Plan is to make Complex(cyc) return Complex{<:QuadraticRing} rather than a float type.
Base.float(cyc::CyclotomicRing{4}) = Complex(cyc)

Base.complex(cyc::CyclotomicRing) = Complex(cyc)

#Base.Complex(cyc::CyclotomicRing{4}) = AbstractFloat(cyc)

# @inline function Base.AbstractFloat(c::CyclotomicRing{D}) where {D}
#     _convert_cyc(c, D, float)
# end

# @inline function Base.Complex(c::CyclotomicRing{D}) where {D}
#     _convert_cyc(c, D, complex)
# end

#function Base.Complex{Tc}(cyc::CyclotomicRing{D}) where {D, Tc}

function Base.Complex(cyc::CyclotomicRing{D, T}) where {D, T}
    Tc = float(T)
    _convert_cyc(cyc, D, Complex{Tc})
end


function Base.Complex{Tc}(cyc::CyclotomicRing{D}) where {D, Tc}
    _convert_cyc(cyc, D, Complex{Tc})
end

function (::Type{T})(cyc::CyclotomicRing{D}) where {T<:Real, D}
    (f, rest...) = coeffs(cyc)
    all(iszero, rest) || throw(ArgumentError(lazy"Inexact error converting $cyc to $T."))
    T(f)
end

Base.float(cyc::CyclotomicRing) = complex(cyc)


# Faster than using any kind of iterative or reduce scheme
# This is only slightly slower 2.7 vs 3.5 ns than the explicit
# expression above for D = 4
@inline function _convert_cyc(c, D, func)
    _add_term(func, 0, D, coeffs(c)...)
end

@inline function _add_term(func, i, D, c)
    func(RootOne{2*D}(i)) * func(c)
end

@inline function _add_term(func, i, D, c, _coeffs...)
    func(RootOne{2*D}(i)) * func(c) + _add_term(func, i+1, D, _coeffs...)
end

# Slightly faster than the general method: 3.6ns vs 2.6ns
function Base.Complex{Tc}(c::CyclotomicRing{4}) where Tc
    (a, b, c, d) = c.coeffs
    r8 = RootOne{8}
    T = Complex{Tc}
    T(a) + T(b) * T(r8(1)) + T(c) * T(r8(2)) +
        T(d) * T(r8(3))
end

# function (::Type{Complex})(c::CyclotomicRing)
#     convert(Complex, c)
# end

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
    CyclotomicRing{4, T}(coeffs)
end

end # module CyclotomicRings
