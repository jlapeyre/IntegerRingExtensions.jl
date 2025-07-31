module CyclotomicRings

import LinearAlgebra
import Base: convert, zero, iszero, one, isone, promote_rule
import ..Common: canonical, imaginary, sqrt_imaginary, one_over_root_two
import ..Utils: superscript
import ..RootOnes: RootOne8

import ..DyadicFractions: DyadicFraction
import ..QuadraticRings: root2conj, Droot2

using ILog2

export CyclotomicRing, Zomega, Domega

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

Domega(a, b, c, d) = Domega{Int}(a,b,c,d)
Domega{T}(a, b, c, d) where {T <: Integer} = _mkomega(DyadicFraction{T,Int}, a,b,c,d)

Zomega(a, b, c, d) = Zomega{Int}(a,b,c,d)
Zomega{T}(a, b, c, d) where {T <: Integer} = _mkomega(T, a,b,c,d)

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

#function Base.:*(c1::CyclotomicRing{4, T1}, c2::CyclotomicRing{4, T2}) where {T1, T2}
function Base.:*(c1::CyclotomicRing{4}, c2::CyclotomicRing{4})
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
    (ai, bi, ci, di) = (c1, d1, -a1, -b1)

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
            (b,c,d,-a)
        elseif k == 2
            (c,d,-a,-b)
        elseif k == 3
            (d,-a,-b,-c)
        elseif k == 4
            (-a,-b,-c,-d)
        elseif k == 5
            (-b,-c,-d,a)
        elseif k == 6
            (-c,-d,a,b)
        elseif k == 7
            (-d,a,b,c)
        end
    return typeof(cyc)(coeffs)
end

Base.:*(c::CyclotomicRing{4}, r::RootOne8) = r * c

@inline Base.:(==)(c1::CyclotomicRing, c2::CyclotomicRing) = c1.coeffs == c2.coeffs

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

function Base.isone(c::CyclotomicRing{4, CoeffT}) where {CoeffT}
    z = zero(CoeffT)
    o = one(CoeffT)
    (a, b, c, d) = c.coeffs
    return a == z && b == z && c == z && d == o
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

function CyclotomicRing{4, T}(r::RootOne8) where {T}
    val = r.k == 4 ? -1 : sign(4 - r.k)
    pos = mod(r.k, 4)
    coeffs =
        if pos == 0
            (0,0,0,val)
        elseif pos == 1
            (0,0,val,0)
        elseif pos == 2
            (0,val,0,0)
        else
            (val,0,0,0)
        end
    CyclotomicRing{4, T}(coeffs)
end

function Base.big(c::CyclotomicRing)
    Complex{BigFloat}(c)
end

function Base.convert(::Type{Complex{Tc}}, c::CyclotomicRing) where Tc
    Complex{Tc}(c)
end

end # module CyclotomicRings
