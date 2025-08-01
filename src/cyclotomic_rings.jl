module CyclotomicRings

import LinearAlgebra
import Base: convert, zero, one, promote_rule
import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two, root_two, coeffs,
    mul_root_two, mul_one_over_root_two, mul_half
import ..Utils: superscript, iszero_strong, isone_strong, PRETTY
import ..RootOnes: RootOne8, RootOne

import ..DyadicFractions: DyadicFraction
import ..QuadraticRings: conj_root_two, Droot2

using ILog2

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
const Domega{T} = CyclotomicRing{4, DyadicFraction{T, Int}}
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

julia> coeffs(Droot2(1,DyadicFraction(3,2)))
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
    # ind = length(cyc.coeffs) - n
    # cyc.coeffs[ind]
    cyc.coeffs[n+1]
    # eh checkbounds does not work for tuple
    # checkbounds(cyc.coeffs, n + 1)
    # @inbounds cyc.coeffs[n+1]
end

# Iterate over coefficients in reverse order.
# Might be better to store them in the usual order (insted of reverse)
# and to display them in reverse order to agree with the literature.
# function Base.iterate(cyc::CyclotomicRing, i::Integer=1)
#     n = length(cyc.coeffs)
#     i > n && return nothing
#     return (cyc.coeffs[n-i+1], i + 1)
# end

function Base.iterate(cyc::CyclotomicRing, i::Integer=1)
    iterate(cyc.coeffs, i)
    # n = length(cyc.coeffs)
    # i > n && return nothing
    # return (cyc.coeffs[n-i+1], i + 1)
end

Base.length(cyc::CyclotomicRing) = length(cyc.coeffs)

Base.lastindex(::CyclotomicRing{M}) where M = M - 1

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

function Base.conj(cyc::Domega{T}) where T
    (a, b, c, d) = cyc.coeffs
    # I think promotion not needed.
#   newcoeff = promote(-c, -b, -a, d)
    newcoeff = promote(a, -d, -c, -b)
    Domega{T}(newcoeff)
end

"""
    conj(cyc::CyclotomicRing{4})

Return the complex conjugate of `cyc`.
"""
function Base.conj(cyc::CyclotomicRing{4})
    (a, b, c, d) = cyc.coeffs
#    newcoeff = (-c, -b, -a, d)
    newcoeff = (a, -d, -c, -b)
    typeof(cyc)(newcoeff)
end

Base.adjoint(cyc::CyclotomicRing) = conj(cyc)

"""
    conj_root_two(cyc::Domega)

Maps `a + b√2` to `a - b√2` in all coefficients in `cyc`.

This implements √2-conjugation.
"""
function conj_root_two(cyc::Domega{T}) where T
    (a, b, c, d) = cyc.coeffs
    newcoeff = (a, -b, c, -d)
    Domega{T}(newcoeff)
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

"""
    one_over_root_two(::Type{Domega{T}}) where {T}

Return a value of `Domega{T}` representing the reciprocal of the square root of two.

# Example
```jldoctest
julia> one_over_root_two(Domega{Int})
1/2 ω + -1/2 ω³
```
"""
function one_over_root_two(::Type{CyclotomicRing{4, DyadicFraction{T, Int}}}) where {T}
    DT = DyadicFraction{T, Int}
    z = zero(DT)
    half = DyadicFraction(T(1), 1)
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
function root_two(::Type{CyclotomicRing{4, DyadicFraction{T, Int}}}) where {T}
    DT = DyadicFraction{T, Int}
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
    CyclotomicRing(cs)
end

# function Domega(c1, coeffs...)
#     cs = promote(DyadicFraction(c1), coeffs...)
#     CyclotomicRing(cs)
# end

function Domega(c1, c2, c3, c4)
    cs = promote(DyadicFraction(c1), c2, c3, c4)
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

Domega{T}(a, b, c, d) where {T <: Integer} = _mkomega(DyadicFraction{T,Int}, a,b,c,d)

Zomega(a, b, c, d) = Zomega{Int}(a,b,c,d)
Zomega{T}(a, b, c, d) where {T <: Integer} = _mkomega(T, a,b,c,d)

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
        a1*a2 - c1*c2 - d1*b2 - d2*b1,

        b2*a1 + b1*a2 - c1*d2 - d1*c2,

        c2*a1 + b1*b2 + c1*a2 - d1*d2,

        d2*a1 + b1*c2 + c1*b2 + d1*a2,
        # a2*d1 + c1*b2 + b1*c2 + a1*d2,
        # b2*d1 + c1*c2 + b1*d2 - a1*a2,
        # c2*d1 + c1*d2 - b1*a2 - a1*b2,
        # d1*d2 - b1*b2 - a1*c2 - a2*c1

    )
    CyclotomicRing(coeffs)
end

function Base.:*(c::CyclotomicRing, r::Real)
    cs = c.coeffs .* r
    CyclotomicRing(cs)
end

@inline Base.:*(r::Real, c::CyclotomicRing) = c * r

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

function Base.:*(r::RootOne8, cyc::CyclotomicRing{4})
    k = r.k
    (a, b, c, d) = cyc.coeffs
    coeffs =
        if k == 0
            (a,b,c,d)
        elseif k == 1
            (-d,a,b,c)
        elseif k == 2
            (-c,-d,a,b)
        elseif k == 3
            (-b,-c,-d,a)
        elseif k == 4
            (-a,-b,-c,-d)
        elseif k == 5
            (d,-a,-b,-c)
        elseif k == 6
            (c,d,-a,-b)
        elseif k == 7
            (b,c,d,-a)
        end
        # if k == 0
        #     (a,b,c,d)
        # elseif k == 1
        #     (b,c,d,-a)
        # elseif k == 2
        #     (c,d,-a,-b)
        # elseif k == 3
        #     (d,-a,-b,-c)
        # elseif k == 4
        #     (-a,-b,-c,-d)
        # elseif k == 5
        #     (-b,-c,-d,a)
        # elseif k == 6
        #     (-c,-d,a,b)
        # elseif k == 7
        #     (-d,a,b,c)
        # end
    return typeof(cyc)(coeffs)
end

Base.:*(c::CyclotomicRing{4}, r::RootOne8) = r * c

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

"""
    mul_half(cyc::CyclotomicRing{4})

Multiply `cyc` by the reciprocal of two.
"""
function mul_half(cyc::CyclotomicRing{4}, n::Integer=1)
    CyclotomicRing(map(x -> mul_half(x, n), cyc.coeffs))
end

"""
    mul_one_over_root_two(cyc::CyclotomicRing{4})

Multiply `cyc` by the reciprocal of the square root of two.
"""
function mul_one_over_root_two(cyc::CyclotomicRing{4})
    mul_root_two(mul_half(cyc))
end

@inline Base.:(==)(c1::CyclotomicRing, c2::CyclotomicRing) = c1.coeffs == c2.coeffs

# @inline function Base.:(==)(c1::CyclotomicRing{N}, r::RootOne{M}) where {N, M}
#     isone(c1) && isone(r) && return true
#     return false
# end

function Base.:^(c::CyclotomicRing, n::Integer)
    n == 0 && return one(c)
    n == 1 && return c
    n == 2 && return c * c
    return Base.power_by_squaring(c, n)
end

function Base.one(::Type{CyclotomicRing{D, CoeffT}}) where {D, CoeffT}
    CyclotomicRing{D}(one(CoeffT), ntuple(x -> zero(CoeffT), D-1)...,)
end

Base.one(c::CyclotomicRing) = one(typeof(c))

Base.zero(::Type{CyclotomicRing{D, CoeffT}}) where {D, CoeffT} = CyclotomicRing(ntuple(x -> zero(CoeffT), D)...,)
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
    float(a) + float(b) * float(RootOne8(1)) + float(c) * float(RootOne8(2)) +
        float(d) * float(RootOne8(3))
end

@inline function Base.AbstractFloat(c::CyclotomicRing{D}) where {D}
    _convert_cyc(c, D, float)
end

@inline function Base.complex(c::CyclotomicRing{D}) where {D}
    _convert_cyc(c, D, complex)
end

function Base.Complex{Tc}(c::CyclotomicRing{D}) where {D, Tc}
    _convert_cyc(c, D, Complex{Tc})
end

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
    T = Complex{Tc}
    T(a) + T(b) * T(RootOne8(1)) + T(c) * T(RootOne8(2)) +
        T(d) * T(RootOne8(3))
end

function CyclotomicRing{4, T}(r::RootOne8) where {T}
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
        # if pos == 0
        #     (0,0,0,val)
        # elseif pos == 1
        #     (0,0,val,0)
        # elseif pos == 2
        #     (0,val,0,0)
        # else
        #     (val,0,0,0)
        # end
    CyclotomicRing{4, T}(coeffs)
end

function Base.big(c::CyclotomicRing)
    Complex{BigFloat}(c)
end

function Base.convert(::Type{Complex{Tc}}, c::CyclotomicRing) where Tc
    Complex{Tc}(c)
end

end # module CyclotomicRings
