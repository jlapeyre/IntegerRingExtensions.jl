module RingMatrices

import ..CyclotomicRings: DOmega, ZOmega, least_denominator_exponent
import ..Matrices2x2: Matrix2x2, AbstractMatrix2x2, elements, ScaleMatrix2x2
import ..RootOnes: RootOne
import IsApprox: isinvolution, Approx, AbstractApprox, Equal

import ..Singletons: InvTwo, InvTwoT,
    RootTwo, RootTwoT, Two, TwoT,
    InvRootTwo, InvRootTwoT,
    Imag, ImagT,
    RootImag, RootImagT,
    Pow

# The modules for rings, groups, and Matrices are independent.
# This module mixes them.

function Base.:*(r::RootOne{8}, m::Matrix2x2{<:DOmega})
    map(x -> r * x, m)
end

Base.:*(m::Matrix2x2{<:DOmega}, r::RootOne{8}) = r * m

function compute_phase_factor(m1::Matrix2x2{<:DOmega}, m2::Matrix2x2{<:DOmega})
    for i in 0:7
        RootOne{8}(i) * m1 == m2 && return i
    end
    return -1
end

function least_denominator_exponent(m::AbstractMatrix2x2{<:DOmega})
    maximum(least_denominator_exponent, elements(m))
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

end #module RingMatrices
