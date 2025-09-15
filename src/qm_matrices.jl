module QMMatrices

import LinearAlgebra
import LinearAlgebra: eigvals, svdvals, opnorm, tr, det, diag, diagm, eigvecs, norm, normalize
import IsApprox:  AbstractApprox, Equal, Approx, isposdef, ispossemidef, isidempotent # isunitary, isinvolution

import ..Matrices2x2: AbstractMatrix2x2, Matrix2x2, AbstractNormalNxN, elements

const I2x2 = Matrix2x2(1, 0, 0, 1)
const X = Matrix2x2(0, 1, 1, 0)
const Y = Matrix2x2(complex(0), 1im, -1im, complex(0))
const Z = Matrix2x2(1, 0, 0, -1)

random_density_matrix2x2(;pure=false) = random_density_matrix2x2(Float64;pure=pure)
function random_density_matrix2x2(::Type{T} ;pure=false) where {T <: AbstractFloat}
    (px, py, pz) = rand(T, 3)
    s = sqrt(px*px + py*py + pz*pz)
    pure || (s = s / rand(T))
    (px, py, pz) = (px, py, pz) ./ s
    return 0.5 * (I2x2 + px * X + py * Y + pz * Z)
end

abstract type AbstractDensityMatrixNxN{T, N} <: AbstractNormalNxN{Complex{T}, N} end

const AbstractDensityMatrix2x2{T} = AbstractDensityMatrixNxN{T, 2} where {T}

# We can try to not require AbstractFloat, so that we can use symbolic matrices
struct DensityMatrix2x2{T} <: AbstractDensityMatrix2x2{T}
    px::T
    py::T
    pz::T
end

struct PureDensityMatrix2x2{T<:AbstractFloat} <: AbstractDensityMatrix2x2{T}
    px::T
    py::T
end

isdensitymatrix(::AbstractDensityMatrixNxN, a::AbstractApprox=Equal()) = true
@inline isdensitymatrix(m::AbstractMatrix2x2) = isdensitymatrix(m, Equal())
@inline function isdensitymatrix(m::AbstractMatrix2x2, approx::AbstractApprox)
    ispossemidef(m, approx) || return false
    isapprox(tr(m), 1, approx)
end

polarization_norm(rho::PureDensityMatrix2x2) = 1

function polarization(rho::PureDensityMatrix2x2)
    (px, py) = (rho.px, rho.py)
    pz = sqrt(1 - px*px - py*py)
    (px, py, pz)
end

# Need a sampler, and interface, etc.
# function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{DensityMatrix2x2})
#     DensityMatrix2x2(random_density_matrix2x2(rng))
# end

"""
    polarization(rho::DensityMatrix2x2)

Return the polarization vector of `rho`.
"""
polarization(rho::DensityMatrix2x2) = (rho.px, rho.py, rho.pz)

"""
    polarization_norm(rho::DensityMatrix2x2)

Return the square of the norm of the polarization vector of `rho`.
"""
function polarization_norm(rho::DensityMatrix2x2)
    (x, y, z) = polarization(rho)
    x*x + y*y + z*z
end

function DensityMatrix2x2(m::AbstractMatrix2x2; check::Bool=false, approx::AbstractApprox = Equal())
    if check || !isa(approx, Equal)
        isdensitymatrix(m, approx) || throw(ArgumentError(lazy"Argument is not a density matrix"))
    end
    DensityMatrix2x2(real(tr(X * m)), real(tr(Y * m)), real(tr(Z * m)))
end

function PureDensityMatrix2x2(m::AbstractMatrix2x2)
   PureDensityMatrix2x2(real(tr(X * m)), real(tr(Y * m)))
end

function elements(m::PureDensityMatrix2x2)
    elements(DensityMatrix2x2(m))
end

function DensityMatrix2x2(r::PureDensityMatrix2x2)
    DensityMatrix2x2(polarization(r)...)
end

function Matrix2x2(r::PureDensityMatrix2x2)
    Matrix2x2(DensityMatrix2x2(r))
end

@inline _m00(r::DensityMatrix2x2) = (1+r.pz)/2
@inline _m01(r::DensityMatrix2x2) = (r.px - r.py * im)/2
@inline _m10(r::DensityMatrix2x2) = (r.px + r.py * im)/2
@inline _m11(r::DensityMatrix2x2) = (1-r.pz)/2

@inline function Matrix2x2(rho::AbstractDensityMatrix2x2{T}) where {T <: AbstractFloat}
    Matrix2x2{Complex{T}}(rho)
end

@inline function Matrix2x2{Complex{T}}(rho::DensityMatrix2x2) where {T<:AbstractFloat}
#    (px, py, pz) = polarization(rho)
    CT = Complex{T}
    Matrix2x2(CT(_m00(rho)), CT(_m10(rho)), CT(_m01(rho)), CT(_m11(rho)))
end

function Base.getindex(rho::DensityMatrix2x2, i::Integer)
    i == 1 && return complex(_m00(rho))
    i == 2 && return _m10(rho)
    i == 3 && return _m01(rho)
    i == 4 && return complex(_m11(rho))
    throw(BoundsError(rho, i))
end

tr(rho::AbstractDensityMatrix2x2) = 1

# All unitary matrices are normal. All Hermitian matrices are normal
# Not all Hermitian matrices are unitary. So matrices are not categorized
# into a tree by these features.
# So it would be useful to implement traits for these features.
LinearAlgebra.ishermitian(rho::AbstractDensityMatrix2x2) = true
LinearAlgebra.ishermitian(rho::AbstractDensityMatrix2x2, ::AbstractApprox) = true
LinearAlgebra.ishermitian(rho::AbstractDensityMatrix2x2, ::Approx) = true
Base.adjoint(rho::AbstractDensityMatrix2x2) = rho

diag(rho::DensityMatrix2x2) = Vector2(real(1+rho.pz)/2, real(1-rho.pz)/2)

ispure(::PureDensityMatrix2x2) = true
ispure(::PureDensityMatrix2x2, ::AbstractApprox) = true
ispure(rho::DensityMatrix2x2) = ispure(rho, Equal())
function ispure(rho::DensityMatrix2x2, approx::AbstractApprox)
    isapprox(1 + polarization_norm(rho), 2, approx)
end

function ispure(m::AbstractMatrix, approx::AbstractApprox=Equal())
    isdensitymatrix(m, approx)
    a = @inbounds m[1,1]
    b = @inbounds m[1,2]
    isapprox(a*a + (1 - a)*(1 - a) + 2*abs2(b), 1, approx)
end

det(::PureDensityMatrix2x2) = 0
function det(rho::DensityMatrix2x2)
    1 - polarization_norm(rho)
end

end # module QMMatrices
