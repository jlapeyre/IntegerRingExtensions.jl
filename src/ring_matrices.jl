module RingMatrices

import LinearAlgebra
import IsApprox: isinvolution, Approx, AbstractApprox, Equal
import ..CyclotomicRings: DOmega, ZOmega, least_denominator_exponent, CyclotomicRing, div_coefficients
import ..Matrices2x2: MatrixNxN, Matrix2x2, AbstractMatrix2x2, AbstractMatrixNxN, elements, ScaleMatrix2x2, scalematrix
import ..RootOnes: RootOne, Omega
import ..Common: canonical, coeffs
import ..Utils: lobit

import ..Singletons: InvTwo, InvTwoT,
    RootTwo, RootTwoT, Two, TwoT,
    InvRootTwo, InvRootTwoT,
    Imag, ImagT,
    RootImag, RootImagT,
    Pow, SingleNum

# The modules for rings and groups on one hand, and `MatrixNxN` on the other, are independent.
# This module mixes them.

function Base.:*(r::RootOne{8}, m::MatrixNxN{<:CyclotomicRing{4}})
    map(x -> r * x, m)
end

Base.:*(m::AbstractMatrixNxN, r::RootOne{8}) = throw(MethodError(Base.:*, (m, r)))

# Organize better. These should only be more generic methods, in say, matrices.jla
Base.:*(m::AbstractMatrixNxN{<:ZOmega}, r::RootOne{8}) = r * m
Base.:*(m::AbstractMatrixNxN{<:DOmega}, r::RootOne{8}) = r * m

Base.:*(n::SingleNum, m::AbstractMatrixNxN) = map(x -> n * x, m)
Base.:*(m::AbstractMatrixNxN, n::SingleNum) = map(x -> n * x, m)

"""
    compute_phase_factor(m1::Matrix2x2{<:DOmega}, m2::Matrix2x2{<:DOmega})::Int

Find `k` in `0:7` such that `Omega(k) * m1 == m2`. If no such `k` exists, return `-1`.
"""
function compute_phase_factor(m1::Matrix2x2{<:DOmega}, m2::Matrix2x2{<:DOmega})
    for i in 0:7
        RootOne{8}(i) * m1 == m2 && return i
    end
    return -1
end

"""
    least_denominator_exponent(m::AbstractMatrix2x2{<:DOmega})

The minimum `k`, such that `√2ᵏ * m[i,j] ∈ ℤ[ω]` for all `i,j`.
"""
function least_denominator_exponent(m::AbstractMatrix2x2{<:DOmega})
    maximum(least_denominator_exponent, elements(m))
end

# This is, of course, very fast. Because the lde is stored in a field.
function least_denominator_exponent(m::ScaleMatrix2x2{<: DOmega, <: Matrix2x2{<:ZOmega}, <: Pow{InvRootTwoT}})
    m.s.n
end

isinvolution(m::AbstractMatrix{<:DOmega}) = isinvolution(m, Equal())
function isinvolution(m::AbstractMatrix{<:DOmega}, ::Equal)
    isone(m * m)
end

function isinvolution(m::AbstractMatrix{<:DOmega}, approx::Approx)
    isone(m * m, approx)
end

"""
    scalematrix(m::AbstractMatrix{<:DOmega})::ScaleMatrix2x2

Convert `m` to a `ScaleMatrix2x2` containing a matrix of type `Matrix2x2{<:ZOmega}` and a global power of `1/√2`.
"""
function scalematrix(m::AbstractMatrix{<:DOmega})
    k = least_denominator_exponent(m)
    m = RootTwo^k * m
    T = typeof(ZOmega(m[1]))
    mz = Matrix2x2{T}(m)
    ScaleMatrix2x2(mz, InvRootTwo^k)
end

function scalematrix(m::AbstractMatrix{<:ZOmega})
    ScaleMatrix2x2(m, InvRootTwo^0)
end

# Hmm this does not depend on rings at all
# function canonical(sm::ScaleMatrix2x2)
#     ScaleMatrix2x2(canonical(sm.m), sm.s)
# end

# Keep getting DOmega as the matrix inside scale matrix
function canonical(sm::ScaleMatrix2x2{<:DOmega, <: Any, <: Any})
    throw(MethodError(canonical, (sm,)))
end

# This seems to work, but I did not test much, because it is a PITA
# to integrate into compose
function canonical(sm::ScaleMatrix2x2{<:DOmega, Matrix2x2{V}, Pow{T}}) where {T, V <: ZOmega}
    # Number of inverse factors of two in the coefficient.
    num_twos_coeff = div(sm.s.n, 2) # Assume this is pow of sqrt(2)
    iszero(num_twos_coeff) && return sm
    m = canonical(sm.m)
    # Smallest power of two in factorization over coefficients and elements of matrix
    minpowtwo = minimum(xx -> minimum(lobit, coeffs(xx)), m)
    if minpowtwo >= num_twos_coeff
        newm = map(x -> div_coefficients(x, 2^num_twos_coeff), m)
        new_s = sm.s * Two^num_twos_coeff
        return ScaleMatrix2x2(newm, new_s)
    end
    new_s = Two^minpowtwo * sm.s
    newm = map(x -> div_coefficients(x, 2^minpowtwo), m)
    @assert !isa(newm, Matrix2x2{<:DOmega})
    return ScaleMatrix2x2(newm, new_s)
end

function Base.:*(r::RootOne{8}, sm::ScaleMatrix2x2{<:DOmega})
    #    typeof(sm)(r * sm.m, sm.s) # TODO: write method for this.
    ScaleMatrix2x2(r * sm.m, sm.s)
end

# Convert m to SU2 by multiplying second column by a power of omega.
# We could represent unitary m by this SU2 and the power.
function getsu2(m::Matrix2x2{<:CyclotomicRing{4}})
    ph = inv(Omega(LinearAlgebra.det(m)))
    (a, b, c, d) = elements(m)
    Matrix2x2(a, b, ph * c, ph * d)
end

end #module RingMatrices
