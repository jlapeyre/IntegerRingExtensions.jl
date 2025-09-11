module RingMatrices

import ..CyclotomicRings: DOmega, ZOmega, least_denominator_exponent
import ..Matrices2x2: MatrixNxN, Matrix2x2, AbstractMatrix2x2, AbstractMatrixNxN, elements, ScaleMatrix2x2
import ..RootOnes: RootOne
import ..Common: canonical
import IsApprox: isinvolution, Approx, AbstractApprox, Equal


import ..Singletons: InvTwo, InvTwoT,
    RootTwo, RootTwoT, Two, TwoT,
    InvRootTwo, InvRootTwoT,
    Imag, ImagT,
    RootImag, RootImagT,
    Pow, SingleNum

# The modules for rings and groups on one hand, and `MatrixNxN` on the other, are independent.
# This module mixes them.

function Base.:*(r::RootOne{8}, m::MatrixNxN{<:DOmega})
    map(x -> r * x, m)
end

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
#    @show typeof(mz)
    ScaleMatrix2x2(mz, InvRootTwo^k)
end

function canonical(sm::ScaleMatrix2x2)
    ScaleMatrix2x2(canonical(sm.m), sm.s)
end

function Base.:*(r::RootOne{8}, sm::ScaleMatrix2x2{<:DOmega})
    #    typeof(sm)(r * sm.m, sm.s) # TODO: write method for this.
    ScaleMatrix2x2(r * sm.m, sm.s)
end


end #module RingMatrices
