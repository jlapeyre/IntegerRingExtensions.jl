module CyclotomicRings

import LinearAlgebra
import Base: convert, zero, iszero, one, isone, promote_rule
import ..IntegerExtensions: imaginary, sqrt_imaginary, one_over_root_two
import ..Common: canonical
import ..Utils: superscript
import ..RootOnes: RootOne8

import ..DyadicFractions: DyadicFraction
import ..QuadraticRings: root2conj

using ILog2

export Domega, CyclotomicRing

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

function Base.:*(c::CyclotomicRing, z::Complex)
    rs = real(z) .* c.coeffs
    (a, b, c, d) = imag(z) .* c.coeffs
    cs = rs .+ (c, d, -a, -b)
    CyclotomicRing(cs)
end

Base.:*(z::Complex, c::CyclotomicRing) = c * z


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

end # module CyclotomicRings
