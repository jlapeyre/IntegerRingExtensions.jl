module Matrices2x2
#@stable module Matrices2x2

# using DispatchDoctor: @stable

import Random
import LinearAlgebra: eigvals, svdvals, opnorm, tr, det, diag, diagm, eigvecs, eigen, norm, normalize,
    dot, schur, Schur, isdiag
import LinearAlgebra: LinearAlgebra, isposdef, ishermitian, hermitianpart, hermitian, issymmetric
import IsApprox: isunitary, isinvolution, AbstractApprox, Equal, Approx, ispossemidef,
    isnormal

import ..Common: canonical, isunit, isimag
import ..Utils: PRETTY, cpad, _show_with_fieldnames, _power_by_squaring
import ..Angles: radtodar, Dar, Ang, intdiv, random_angle

export Matrix2x2, AbstractMatrix2x2, Matrix4x4, AbstractMatrix4x4, antihermitianpart,
    isantihermitian

abstract type AbstractMatrixNxN{T, N} <: AbstractMatrix{T} end
abstract type AbstractNormalNxN{T, N} <: AbstractMatrixNxN{T, N} end
abstract type AbstractUnitaryNxN{T, N} <: AbstractNormalNxN{T, N} end

const AbstractMatrix2x2{T} = AbstractMatrixNxN{T, 2} where T
const AbstractNormal2x2{T} = AbstractNormalNxN{T, 2} where {T}
const AbstractUnitary2x2{T} = AbstractUnitaryNxN{T, 2} where {T}

const AbstractMatrix4x4{T} = AbstractMatrixNxN{T, 4} where T
const AbstractNormal4x4{T} = AbstractNormalNxN{T, 4} where {T}
const AbstractUnitary4x4{T} = AbstractUnitaryNxN{T, 4} where {T}

abstract type AbstractSU2{T} <: AbstractUnitary2x2{T} end
abstract type AbstractVector2{T} <: AbstractVector{T} end

"""
    Matrix2x2{T} <: AbstractMatrix{T}

Stack allocated, immutable, 2 x 2 matrix.

If `T` is `isbits`, then `Matrix2x2{T}` is `isbits`.

We implement a few necessary operations, including matrix multiplication,
addition, subtraction, and unary minus.
"""
struct MatrixNxN{T, N, N2} <: AbstractMatrixNxN{T, N}
    function MatrixNxN{T, N}(data::NTuple{N2, T}) where {T, N, N2}
        N2 == N * N ||  throw(ArgumentError(lazy"Inexact error bad tup"))
        new{T, N, N2}(data)
    end

    function MatrixNxN{T, N, N2}(data::NTuple{N2, T}) where {T, N, N2}
        N2 == N * N ||  throw(ArgumentError(lazy"Inexact error bad tup"))
        new{T, N, N2}(data)
    end

    data::NTuple{N2, T}
end

const Matrix2x2{T} = MatrixNxN{T, 2, 4} where {T}
const Matrix4x4{T} = MatrixNxN{T, 4, 16} where {T}

function Matrix2x2(tup::NTuple{4, T}) where T
    Matrix2x2{T}(tup)
end

# Thank you Elrod
# https://discourse.julialang.org/t/is-this-test-detect-unbound-args-result-valid-or-a-bug/96987
function MatrixNxN(tup::Tuple{T, Vararg{T, N2mOne}}) where {N2mOne, T}
    N = isqrt(N2mOne+1)
    MatrixNxN{T, N, N2mOne+1}(tup)
end

struct Vector2{T} <: AbstractVector2{T}
    data::NTuple{2, T}
end

function columns(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    (Vector2(a, b), Vector2(c, d))
end

##
## Constructors
##

"""
    Matrix2x2(a, b, c, d)

Return a `Matrix2x2` with the given elements.

`Matrix2x2` is column major, as is the type `Matrix`. The columns are
`[a, b]` and `[c, d]`.
"""
Matrix2x2(a, b, c, d) = Matrix2x2(promote(a, b, c, d))
Matrix2x2{T}(a, b, c, d) where {T} = Matrix2x2(T(a), T(b), T(c), T(d))

function Matrix2x2(m::Matrix)
    MatrixNxN{eltype(m), 2}(m)
end

function Matrix4x4(m::Matrix)
    MatrixNxN{eltype(m), 4}(m)
end

function MatrixNxN{T, N}(m::Matrix) where {T, N}
    size(m) == (N, N) || error("Wrong size")
    N2 = N * N
    # Use NTuple{N2, eltype(m)} to avoid really bad perf.
    MatrixNxN{T, N, N2}(NTuple{N2, eltype(m)}(m))
end

function MatrixNxN{T, N, N2}(m::Matrix) where {T, N, N2}
    size(m) == (N, N) || error("Wrong size")
    # Use NTuple{N2, eltype(m)} to avoid really bad perf.
    MatrixNxN{T, N, N2}(NTuple{N2, eltype(m)}(m))
end

Base.Matrix(m::AbstractMatrix2x2) = reshape([elements(m)...,], (2,2))

function Matrix2x2{T}(m::AbstractMatrix2x2) where {T}
    Matrix2x2(map(x -> T(x), elements(m))...,)
end

Matrix2x2(m::Matrix2x2) = m

function Base.map(f, m::AbstractMatrixNxN)
    MatrixNxN(map(f, elements(m)))
end

elements(m::MatrixNxN) = m.data
elements(m::AbstractMatrixNxN{T,N}) where {T, N} = elements(MatrixNxN{T, N, N*N}(m))

##
## Conversion and related construction
##

Base.convert(::Type{Matrix2x2{T}}, m::Matrix2x2) where {T} = Matrix2x2{T}(m)
Base.AbstractFloat(m::Matrix2x2) = map(float, m)

import Base: real, imag, big, complex, float


"""
    canonical(m::Matrix2x2)

Return a new matrix by calling `canonical` element-wise on `m`.

This may reduce or canonicalize elements that implement `canonical`.
"""
canonical

##
## Display
##

function _showstr(obj)
    b = IOBuffer()
    show(IOContext(b, :compact=>true), PRETTY(), obj)
    return String(take!(b))
end

function Base.show(io::IO, ::PRETTY, m::AbstractMatrix2x2)
    summary(io, m)
    println(io, ":")
    _show_matrix2x2(io, m)
end

function _show_matrix2x2(io::IO, m::AbstractMatrix2x2)
    (as, bs, cs, ds)  = map(_showstr, elements(m))
    (al, bl, cl, dl) = map(length, (as, bs, cs, ds))
    w1 = max(al, bl)
    w2 = max(cl, dl)
    spc = "  "
    print(io, cpad(as, w1), spc)
    println(io, cpad(cs, w2))
    print(io, cpad(bs, w1), spc)
    print(io, cpad(ds, w2))
end

function Base.show(io::IO, m::Matrix2x2)
    print(io, typeof(m), "(")
    i = 0
    for x in elements(m)
        i += 1
        if i > 1
            print(io, ", ")
        end
        print(io, x)
    end
    print(io, ")")
end

function Base.show(io::IO, ::PRETTY, v::AbstractVector2)
    summary(io, v)
    println(io)
    (v1, v2) = elements(v)
    show(io, PRETTY(), v1)
    println(io)
    show(io, PRETTY(), v2)
end

##
## Required, and standard, `Base` properties
##

# Annoying to do this piecemeal!!
# unsustainable
function Base.getindex(m::Matrix2x2, ::Colon, i::Integer)
    if i == 1
        return Vector2(m[1,1], m[2,1])
    elseif i == 2
        return Vector2(m[1,2], m[2,2])
    end
    throw(BoundsError(m, i))
end

function Base.getindex(m::Matrix2x2, i::Integer, ::Colon)
    if i == 1
        return Vector2(m[1,1], m[1,2])
    elseif i == 2
        return Vector2(m[2,1], m[2,2])
    end
    throw(BoundsError(m, i))
end

@inline function Base.getindex(m::MatrixNxN, i::Integer)
    @boundscheck checkbounds(m, i)
    return @inbounds m.data[i]
end

Base.getindex(m::AbstractMatrix2x2, i::Integer) = Matrix2x2(m)[i]

function Base.getindex(m::AbstractSU2, i::Integer)
    i == 1 && return unitary_u(m)
    i == 2 && return unitary_t(m)
    i == 3 && return -conj(unitary_t(m))
    return conj(unitary_u(m))
end

function elements(m::AbstractSU2)
    a = unitary_u(m)
    b = unitary_t(m)
    (a, b, -conj(b), conj(a))
end


Base.size(::AbstractMatrixNxN{T, N}) where{T, N} = (N, N)
Base.eltype(::AbstractMatrixNxN{T, N}) where{T, N} = T
Base.IndexStyle(::Type{<:AbstractMatrixNxN}) = IndexLinear()
Base.one(::Type{Matrix2x2{T}}) where T = Matrix2x2(one(T), zero(T), zero(T), one(T))
Base.one(::Matrix2x2{T}) where {T} = one(Matrix2x2{T})
Base.zero(::Type{Matrix2x2{T}}) where T = Matrix2x2(zero(T), zero(T), zero(T), zero(T))
Base.zero(::Matrix2x2{T}) where {T} = zero(Matrix2x2{T})

# Base.isone fallback is efficient for MatrixNxN

# OK, but prbly not needed
# function Base.isone(m::Matrix2x2)
#     (a, b, c, d) = elements(m)
#     isone(a) && isone(d) && iszero(b) && iszero(c)
# end

Base.iszero(m::MatrixNxN) = all(iszero, elements(m))

# This is too liberal in general. i.e. will be wrong
# This is a necessary condition. But there may exist no
# inverse of the same type
function isunit(m::AbstractMatrixNxN)
    !iszero(det(m))
end

##
## IsApprox.isunitary
##

function isunitary(m::MatrixNxN)
    isone(m * m')
end

isunitary(m::MatrixNxN, ::Equal) = isunitary(m)
isunitary(m::AbstractUnitaryNxN) = true
isnormal(m::AbstractUnitaryNxN) = true
Base.iszero(m::AbstractUnitaryNxN) = false

isinvolution(m::AbstractMatrixNxN) = isinvolution(m, Equal())
function isinvolution(m::AbstractMatrixNxN, ::Equal)
    isone(m * m)
end
function isinvolution(m::AbstractMatrixNxN, approx::Approx)
    isone(m * m, approx)
end

function svdvals(::AbstractUnitary2x2{T}) where {T}
    (one(real(T)), one(real(T)))
end

"""
    isunitary(m::Matrix2x2, app::Approx)

Return `true` if `m` is approximately unitary.

This is problematic if `m` is approximately a rotation matrix.
The off-diagonals will have low precision. You might need to pass,
for example, `Approx(rtol=1e-6)` for `app`. But, this tolerance will apply
 to the diagonals as well, which should have more precision in principle.
"""
function isunitary(m::Matrix2x2, app::Approx)
    app == Equal() && return isunitary(m)
    (a, b, c, d) = elements(m)
    isapprox(abs2(a) + abs2(c), one(a); app.kw...) || return false
    isapprox(abs2(b) + abs2(d), one(a); app.kw...) || return false

    isapprox(conj(a) * b, - d * conj(c); app.kw...)
end

##
## Matrix arithmetic
##

@inline function _mul2(m1, m2)
    (a1, b1, c1, d1) = elements(m1)
    (a2, b2, c2, d2) = elements(m2)
    Matrix2x2(a1*a2 + c1*b2, a2*b1 + d1*b2, a1*c2 + c1*d2, b1*c2 + d1*d2)
end

@inline function Base.:*(m::AbstractMatrixNxN{<:Number, N}, x::Number) where {N}
    els = map(z -> x * z, elements(m))
    MatrixNxN{typeof(first(els)), N}(els)
end

@inline function Base.:+(m::AbstractMatrixNxN{<:Number, N}, x::Number) where {N}
    els = map(z -> x + z, elements(m))
    MatrixNxN{typeof(first(els)), N}(els)
end

@inline Base.:*(x::Number, m::AbstractMatrix2x2) = m * x
@inline Base.:+(x::Number, m::AbstractMatrix2x2) = m + x

@inline function _pair_op(op, m1, m2)
    (a1, b1, c1, d1) = elements(m1)
    (a2, b2, c2, d2) = elements(m2)
    Matrix2x2(op(a1, a2), op(b1, b2), op(c1, c2), op(d1,d2))
end

Base.:*(m1::AbstractMatrix2x2, m2::AbstractMatrix2x2) = _mul2(m1, m2)
Base.:+(m1::AbstractMatrix2x2, m2::AbstractMatrix2x2) = _pair_op(+, m1, m2)
Base.:-(m1::Matrix2x2, m2::Matrix2x2) = _pair_op(-, m1, m2)
Base.:-(m::AbstractMatrix2x2) = map(-, m) # Unary minus

Base.:-(a::AbstractMatrix2x2, b::AbstractMatrix2x2) = Matrix2x2(a) - Matrix2x2(b)

# This is only called if `n` is not literal at the call site.
function Base.:^(m::MatrixNxN, n::Integer)
    n == 0 && return one(m)
    n == 1 && return m  # If elements of m are mutable, this is problematic
    n == 2 && return m * m
    return _power_by_squaring(m, n)
end

Base.literal_pow(::typeof(Base.:^), m::MatrixNxN, ::Val{0}) = one(m)
Base.literal_pow(::typeof(Base.:^), m::MatrixNxN, ::Val{1}) = m # If elements of m are mutable, this is problematic
Base.literal_pow(::typeof(Base.:^), m::MatrixNxN, ::Val{2}) = m * m
Base.literal_pow(::typeof(Base.:^), m::MatrixNxN, ::Val{3}) = m * m * m
Base.literal_pow(::typeof(Base.:^), m::MatrixNxN, ::Val{4}) = (m * m) * (m * m)

##
## Linear algebra and related operations
##

function Base.permutedims(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    Matrix2x2(a, c, b, d)
end

Base.adjoint(m::Matrix2x2) = permutedims(map(adjoint, m))
Base.transpose(m::Matrix2x2) = permutedims(map(transpose, m))

@inline function tr(m::AbstractMatrix4x4)
    sum(i -> m[i,i], 1:4)
end

# OK. But prev function is fine
# @inline function tr(m::AbstractMatrix2x2)
#     (a, _b, _c, d) = elements(m)
#     a + d
# end

@inline function det(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    canonical(a * d - c * b)
end

@inline tr(m::AbstractSU2) = 2*real(m[1])

# Ugh. Going through gymnastics to try to get small union splitting.
# It works in LinearAlgebra because `Vector` is returned.
# These `Tuple`s are more complicated.
#
"""
    eigvals(m::Matrix2x2)::NTuple{2}

Return a tuple of the eigenvalues of `m`.
"""
function eigvals(m::AbstractMatrix2x2)
    if LinearAlgebra.isdiag(m)
        (v1, v2) = (m[1], m[4])
        return real(v1) > real(v2) ? (v1, v2) : (v2, v1)
    end
    (a, b, c, d) = elements(m)
    amd = a - d
    apd = a + d
    d = amd * amd + 4*b*c
    _evals(apd, d)
#    d >= 0 ? _evals(apd, d) : _evals(complex(apd), complex(d))
end

function _eigvals(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    amd = a - d
    apd = a + d
    discr_arg = amd * amd + 4*b*c
    _evals(apd, discr_arg)
end

function eigvals(m::AbstractMatrix2x2{<:Complex})
    _eigvals(m)
end

function _evals(apd::Real, darg::Real)
    if darg < 0
        return _evals(complex(apd), complex(darg))
    end
    disc = sqrt(darg)
    ((apd + disc)/2, (apd - disc)/2)
end

function _evals(apd::Complex, darg::Complex)
    disc = _csqrt(darg)
    ((apd + disc)/2, (apd - disc)/2)
end

function _csqrt(x::Real)
    x < 0 ? sqrt(complex(x)) : sqrt(x)
end

function _csqrt(x::Complex)
    sqrt(x)
end


"""
    hermitianpart(m::Matrix2x2)::Matrix2x2

Return `(m + m')/2`
"""
function hermitianpart(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    b2 = (b + conj(c)) / 2
    Matrix2x2(real(a), b2, conj(b2), real(d))
end

"""
    antihermitianpart(m::Matrix2x2)::Matrix2x2

Return `(m - m')/2`
"""
function antihermitianpart(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    b2 = (b - conj(c)) / 2
    Matrix2x2(im * imag(a), b2, -conj(b2), im * imag(d))
end

"""
    hermitian(m::Matrix2x2)::Matrix2x2

Return an instance of `m` converted to a Hermitian matrix.

Because `Matrix2x2` is small and immutable, a new instance is returned rather
than a view. The imaginary parts of `m[1,1]` and `m[2,2]` are set to zero.
`m[2,1]` is set to `conj(m[1,2])`.
"""
function hermitian(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    Matrix2x2(real(a), conj(c), c, real(d))
end

"""
    eigvals_hermitian(m::Matrix2x2)::NTuple{2}

Return a tuple of the eigenvalues of `hermitianpart(m)`

Return type is `<: Real`. If `m` is Hermitian, then `hermitianpart(m)`
is essientially the identity.
"""
function eigvals_hermitian(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(hermitianpart(Matrix2x2(m)))
    amd = real(a - d)
    apd = real(a + d)
    discr = sqrt(amd * amd + 4*abs2(b))
    ((apd + discr)/2, (apd - discr)/2)
end

function eigvals_hermitian_traceone(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(hermitianpart(m))
    amd = real(2*a - 1)
    apd = one(real(a))
    discr = sqrt(amd * amd + 4*abs2(b))
    ((apd + discr)/2, (apd - discr)/2)
end

function LinearAlgebra.norm(v::AbstractVector2)
    (x, y) = elements(v)
    sqrt(abs2(x) + abs2(y))
end

# This fails if the element type changes.
# Typical issue. Workaround is... a bit of work
function LinearAlgebra.normalize(v::AbstractVector2{T}) where {T}
    (x, y) = elements(v)
    n = norm(v)
    (x1, y1) = (x/n, y/n)
    typeof(v)(x/n, y/n)
end

function LinearAlgebra.normalize(v::Vector2{T}) where {T}
    (x, y) = elements(v)
    n = norm(v)
    (x1, y1) = (x/n, y/n)
    Vector2{typeof(x1)}(x1, y1)
end

# Several times faster than generic fallback

"""
    ispossemidef(m::AbstractMatrix2x2, approx::AbstractApprox=Equal())

Return `true` if `m` is positive semi-definite.
"""
@inline ispossemidef(m::AbstractMatrix2x2, approx::AbstractApprox) = _isposdef(m, approx, >=)
@inline ispossemidef(m::AbstractMatrix2x2) = ispossemidef(m, Equal())
@inline ispossemidef(m::AbstractMatrix2x2, approx::Equal) = _isposdef(m, approx, >=)

"""
    isposdef(m::AbstractMatrix2x2, approx::AbstractApprox=Equal())

Return `true` if `m` is positive definite.
"""
@inline isposdef(m::AbstractMatrix2x2, approx::AbstractApprox) = _isposdef(m, approx, >)
@inline isposdef(m::AbstractMatrix2x2) = isposdef(m, Equal())
@inline isposdef(m::AbstractMatrix2x2, approx::Equal) = _isposdef(m, approx, >)

# Not straightforward to approx compare to zero
@inline function _isposdef(m, approx, cmpf::F) where {F}
    if ishermitian(m) # take fast path if possible
        (ev1, ev2) = eigvals_hermitian(m)
        return cmpf(ev1, 0) && cmpf(ev2, 0)
    end
    ishermitian(m, approx) || return false
    (ev1, ev2) = eigvals(m)
    return cmpf(real(ev1), 0) && cmpf(real(ev2), 0)
end

ishermitian(m::AbstractMatrix2x2, approx::AbstractApprox=Equal()) = _ishermitian(m, approx)
ishermitian(m::AbstractMatrix2x2, approx::Approx) = _ishermitian(m, approx)

# We should probably be using Frobenius norm here rather than
# testing elements individually.
# Generic fallback is not less efficient, if we are not using Approx.
# But, with approx, this specialized method is more performant.
function _ishermitian(m::AbstractMatrix2x2, approx::AbstractApprox)
    (a,b,c,d) = elements(m)
    (isreal(a, approx) && isreal(d, approx)) || return false
    isapprox(c, adjoint(b), approx)
end

isantihermitian(m::AbstractMatrix2x2, approx::AbstractApprox=Equal()) = _isantihermitian(m, approx)
isantihermitian(m::AbstractMatrix2x2, approx::Approx) = _isantihermitian(m, approx)
function _isantihermitian(m::AbstractMatrix2x2, approx::AbstractApprox)
    (a,b,c,d) = elements(m)
    (isimag(a, approx) && isimag(d, approx)) || return false
    isapprox(c, -adjoint(b), approx)
end

issymmetric(m::AbstractMatrix2x2, approx::AbstractApprox=Equal()) = _issymmetric(m, approx)
issymmetric(m::AbstractMatrix2x2, approx::Approx) = _issymmetric(m, approx)

function _issymmetric(m::AbstractMatrix2x2, approx::AbstractApprox)
    (b,c) = (m[2], m[3])
    isapprox(c, b, approx)
end

"""
    check_eigen(m, eg::Eigen)

Check that `eigen(m)` is approximately equal to `eg`.
"""
function check_eigen(m, eg::LinearAlgebra.Eigen)
    check_eigen(m, eg.values, eg.vectors)
end

"""
    check_eigen(m, vals, vecs)

Check that `eigen(m)` yields approximately `vals` and `vecs`.
"""
function check_eigen(m, vals, vecs)
    v1 = vecs[:,1]
    v2 = vecs[:,2]
    isapprox(m * v1, vals[1] * v1) &&
        isapprox(m * v2, vals[2] * v2)
end

@inline function _maybe_normalize(v::Vector2)
    iszero(v) ? v : normalize(v)
end

"""
    eigen(m::AbstractMatrix2x2)::Eigen

Return eigenvalues and eigenvectors of `m`.
"""
function eigen(m::AbstractMatrix2x2)
    LinearAlgebra.isdiag(m) && return _eigen_diag(m)
    (v1, v2) = eigvals(m)
    __eigen(m, v1, v2)
end

@inline function __eigen(mi::AbstractMatrix2x2{V}, v1::T, v2::T) where {V, T}
    m = map(typeof(v1), mi)
    (a, b, c, d) = elements(m)
    if (iszero(c) && iszero(v1 - a)) ||
        (iszero(b) && iszero(v2 - d))
        (v1, v2) = (v2, v1)
    end
    _inner_eigen(m, v1, v2)
end

@inline function _inner_eigen(m::AbstractMatrix2x2{T}, v1::T, v2::T) where {T}
    (vec1, vec2) = _get_vecs(m, v1, v2)
    (v11, v21)  = vec1
    (v12, v22)  = vec2
    vtup = _vtup(v11, v21, v12, v22)
    vecs = _get_vec_mat(vtup)
    LinearAlgebra.Eigen(Vector2(v1, v2), vecs)
end

@inline function _vtup(a::T, b::T, c::T, d::T) where {T<:Real}
    (a, b, c, d)
end

@inline function _vtup(a::T, b::T, c::T, d::T) where {T<:Complex}
    (a, b, c, d)
end

@inline function _get_vec_mat(tup::NTuple{4, T}) where {T <: Complex}
    vecs = Matrix2x2{T}(tup)
end

@inline function _get_vec_mat(tup::NTuple{4, T}) where {T <: Real}
    vecs = Matrix2x2{T}(tup)
end

function _get_vecs(m::Matrix2x2{T}, v1::T, v2::T) where {T<:Real}
    (a, b, c, d) = elements(m)
    vec1 = _maybe_normalize(Vector2(c, v1 - a))
    vec2 = _maybe_normalize(Vector2(v2 - d, b))
    if iszero(vec1) && iszero(vec2)
        vec1 = Vector2(one(c), zero(c))
        vec2 = Vector2(zero(c), one(c))
    end
    (vec1, vec2)
end

function _get_vecs(m::Matrix2x2{T}, v1::T, v2::T) where {T<:Complex}
    (a, b, c, d) = elements(m)
    vec1 = _maybe_normalize(Vector2(c, v1 - a))
    vec2 = _maybe_normalize(Vector2(v2 - d, b))
    if iszero(vec1) && iszero(vec2)
        vec1 = Vector2(one(c), zero(c))
        vec2 = Vector2(zero(c), one(c))
    end
    (vec1, vec2)
end

function _eigvals_diag(m)
    (v1, v2) = (m[1], m[4])
    (v1, v2) = real(v1) > real(v2) ? (v1, v2) : (v2, v1)
    Vector2(v1, v2)
end

function _eigen_diag(m)
    evs = _eigvals_diag(m)
    z = zero(evs[1])
    o = one(z)
    return LinearAlgebra.Eigen(evs, Matrix2x2(o, z, z, o))
end

function eigen_hermitian(m::AbstractMatrix2x2)
    (v1, v2) = eigvals_hermitian(m)
    (ai, b, c, di) = elements(m)
    (a, d) = (real(ai), real(di))
    m2 = Matrix2x2(a, b, c, d)
    _inner_eigen(m2, v1, v2)
end

function schur(m::AbstractMatrix2x2)
    (vals, vecs) = eigen(m)
    _schur(vals, vecs, m)
end

function _schur(vals::Vector2{T}, vecs::Matrix2x2{T}, m) where {T}
    v1 = vecs[1]
    v2 = vecs[2]
    U = Matrix2x2(v1, v2, -conj(v2), conj(v1))
    Tm = U' * m * U
    (a, b, c, d) = elements(Tm)
    Tm = Matrix2x2(a, zero(b), c, d)
    return Schur(Tm, U, Vector2((v1, v2)))
end

function _diag_matrix_func(m, func)
    (m1, m4) = (m[1], m[4])
    z = zero(m1)
    Matrix2x2(func(m1), z, z, func(m4))
end

# The wrapper Hermitian can actually expose a Hermitian matrix
# if the input is approx Hermitian.
# But a wrapper Normal... ?
# OTOH, you can assert that it is normal with Normal.
# If it's close enough for you, wrap it in Normal
# If m is normal, compute func(m).
function _normal_matrix_func(m, func)
    isdiag(m) && return _diag_matrix_func(m, func)
    if ishermitian(m)
        ((v1,v2), vecs) = eigen_hermitian(m)
        d = diagm(Vector2(func(complex(v1)), func(complex(v2))))
        result = vecs * d * vecs'
    else
        isnormal(m) || return nothing # AbstractUnitary2x2 takes this path
        ((v1,v2), vecs) = eigen(m)
        D = diagm(Vector2(func(v1), func(v2)))
        result = vecs * D * vecs'
    end
    result
end

function _get_wd(x::T, y::T) where {T}
    ((x+y)/2, (x-y)/2)
end

function secant_line(::typeof(exp), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return exp(w) * sinc(im * d / pi)
end

function secant_line(::typeof(exp), x::T, y::T) where {T <: Real}
    (w, d) = _get_wd(x, y)
    return exp(w) * real(sinc(im * d / pi))
end


function secant_line_full(f::typeof(exp), x::T, y::T) where {T <: Number}
    return (f(x), f(y), secant_line(f, x, y))
end

function secant_line(::typeof(cispi), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return im * pi * cispi(w) * sinc(d)
end

function secant_line_full(::typeof(cispi), x::T, y::T) where {T <: Number}
    return (cispi(x), cispi(y), secant_line(cispi, x, y))
end

function secant_line(::typeof(cis), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return im * cis(w) * sinc(d / pi)
end

function secant_line_full(::typeof(cis), x::T, y::T) where {T <: Number}
    return (cis(x), cis(y), secant_line(cis, x, y))
end

function secant_line(::typeof(sin), x::T, y::T) where {T <: Number}
    (w, d) = _get_wd(x, y)
    return cos(w) * sinc(d / pi)
end

function secant_line_full(::typeof(sin), x::T, y::T) where {T <: Number}
    return (sin(x), sin(y), secant_line(sin, x, y))
end

# secant line of cos
function secant_line(::typeof(cos), x::T, y::T) where {T <: Complex}
    (w, d) = _get_wd(x, y)
    return -sin(w) * sinc(d / pi)
end

function secant_line(::typeof(cos), x::T, y::T) where {T <: Real}
    (w, d) = _get_wd(x, y)
    return -sin(w) * sinc(d / pi)
end

function secant_line_full(::typeof(cos), x::T, y::T) where {T <: Complex}
    return (cos(x), cos(y), secant_line(cos, x, y))
end

function secant_line_full(::typeof(cos), x::T, y::T) where {T <: Real}
    return (cos(x), cos(y), secant_line(cos, x, y))
end

function _approx_diff(f, x, z, h=1e-6)
    w = (x + z)/2
    (f(w+h) - f(w-h)) / (2*h)
end

fastabs(x::Number) = abs(x)
fastabs(z::Complex) = abs(real(z)) + abs(imag(z))

function secant_line_full(f::F, x::T, y::T) where {F, T <: Number}
    fx = f(x)
    fy = f(y)
    dl = x - y
    if fastabs(dl) < sqrt(eps(real(typeof(x))))
        sl = _approx_diff(f, x, y)
    else
        sl = (fx - fy) / (x - y)
    end
    return (fx, fy, sl)
end

# Need this ::F for inference.
# Whether we really need it depends on details of what we return in an unpredicatable way.
function _matrix_func(m::AbstractMatrix2x2{T}, func::F) where {F <: Function} where {T <: Number}
    let normsol = _normal_matrix_func(m, func)
        isnothing(normsol) || return normsol
    end
    (Tm, U, _) = schur(m)
    __matrix_func(U, Tm, func)
end

function __matrix_func(U::T, Tm::T, func::F) where {T, F}
    (x, b, y, z) = elements(Tm)
    (a, d, sline) = secant_line_full(func, x, z)
    c = y * sline
    # @show sline
    # @show a, b, c, d
    tup = _vtup(a, zero(a), c, d)
    md = Matrix2x2(tup)
    return U * md * U'
end

# TODO: Base uses small union splitting to return matrices with Real elements if possible, or else Complex.
# Because we are using matrices backed by Tuple's inference is different. Its a huge PITA
# Small changes can make or break inference.
for func in (:sqrt, :exp, :log, :cos, :sin, :tan, :sinpi, :cospi, :cis, :cispi, :cbrt, :inv,
             :cosh, :sinh, :tanh, :acos, :asin, :atan, :acosh, :asinh, :atanh)
    bfunc = :(Base.$func)
    @eval @inline function $bfunc(m::AbstractMatrix2x2)
        _matrix_func(m, $func)
    end
end

for func in (:cos, :sin, :exp)
    bfunc = :(Base.$func)
    @eval @inline function $bfunc(m::AbstractMatrix2x2{<:Real})
        real(_matrix_func(m, $func))
    end
end

function eigvecs(m::AbstractMatrix2x2)
    (evs, vecs) = eigen(m)
    return vecs
end

# TODO: It might be faster to convert this to Matrix2x2
function eigvals(U::AbstractSU2)
    u = unitary_u(U)
    ru = real(u)
    ri = imag(u)
    abs2t = unitary_abs2t(U)
    rdiscr = sqrt(ri^2 + abs2t)
    discr = Complex(zero(rdiscr), rdiscr)
    (ru + discr, ru - discr)
end

Base.:*(v::AbstractVector2, z::Number) = z * b
function Base.:*(z::Number, v::AbstractVector2)
    (x, y) = elements(v)
    Vector2(z * x, z * y)
end

function LinearAlgebra.dot(A::AbstractVector2, B::AbstractVector2)
    (x, y) = elements(A)
    (v, w) = elements(B)
    x * v' + y * w'
end

"""
    svdvals(m::Matrix2x2)::NTuple{2}

Return a tuple of the singular values of `m` in descending order.
"""
function svdvals(m::AbstractMatrix2x2)
    ma = m * adjoint(m)
    if isone(ma)
        ev = one(real(eltype(m)))
        return (ev, ev)
    end
    (v1, v2) = eigvals(ma)
    (s1, s2) = (sqrt(real(v1)), sqrt(real(v2)))
    s1 > s2 ? (s1, s2) : (s2, s1)
end

function LinearAlgebra.norm(m::AbstractMatrix2x2, n::Real=2)
    n == 2 && return norm2(m)
    n == 1 && return norm1(m)
    n == Inf && return normInf(m)
    throw(ArgumentError(lazy"Unsupported l-norm $n"))
end

# We may need to scale for accuracy, as is done in LinearAlgebra code.
norm2(m::AbstractMatrix2x2) = sqrt(sum(abs2, elements(m)))
norm1(m::AbstractMatrix2x2) = sum(abs, elements(m))
normInf(m::AbstractMatrix2x2) = maximum(abs, elements(m))

function opnorm(m::AbstractMatrix2x2)
    (v1, v2) = svdvals(m)
    max(v1, v2)
end

function opnormdistance(a::AbstractMatrixNxN, b::AbstractMatrixNxN)
    opnorm(a - b)
end

"""
    tracenorm(m::Matrix2x2)

Return the trace norm of `m`.
"""
function tracenorm(m::AbstractMatrixNxN)
    sum(svdvals(m))
 end

function tracedistance(a::AbstractMatrix2x2, b::AbstractMatrix2x2)
    tn = tracenorm(a - b)
    tn / typeof(tn)(2)
end

# This paper calls this the trace distance.
# Authors: Alex Bocharov, Yuri Gurevich, Krysta M. Svore
# http://arxiv.org/abs/1303.1411v1
"""
    GPID(A::Matrix2x2, B::Matrix2x2)

Compute the global phase invariant distance between `A` and `B`.

(See Mukhopadhyay 2021)
"""
function GPID(A::AbstractMatrix2x2, B::AbstractMatrix2x2)
    return sqrt(1 - abs(trace_product(A, B)) / 2)
end

"""
    trace_product(A, B)

Compute tr(AB').

This is the dot product of the matrices as vectors.
"""
function trace_product(A::AbstractMatrix2x2, B::AbstractMatrix2x2)
    (a, b, c, d) = map(conj, elements(A))
    (w, x, y, z) = elements(B)
    return a*w + b*x + c*y + d*z
end

##
## Vector2
##

Vector2(a, b) = Vector2(promote(a, b))
Vector2{T}(a, b) where {T} = Vector2(T(a), T(b))
Base.size(::Vector2) = (2,)
Base.eltype(::Type{Vector2{T}}) where T = T
Base.IndexStyle(::Type{<:Vector2}) = IndexLinear()
elements(v::Vector2) = v.data

@inline function Base.getindex(m::Vector2, i::Integer)
    @boundscheck checkbounds(m, i)
    return @inbounds m.data[i]
end

function Base.map(f, v::Vector2)
    Vector2(map(f, elements(v)))
end

Base.float(v::Vector2) = map(float, v)

function Base.:*(m::Matrix2x2, v::Vector2)
    (a,b,c,d) = elements(m)
    (x, y) = elements(v)
    Vector2(a*x + c*y, b*x + d*y)
end

"""
    diag(m::Matrix2x2)::Vector2

Return the diagonal of `m` as a `Vector2`.
"""
diag(m::Matrix2x2) = Vector2(m[1], m[4])

"""
    diagm(d::Vector2{T})::Matrix2x2 where {T}

Construct a matrix with diagonal given by `d`.
"""
diagm(d::Vector2{T}) where {T} = Matrix2x2(d[1], zero(T), zero(T), d[2])

"""
    random_diagonal_unitary(::Type{T}=Float64)::Matrix2x2 where T

Return a random diagonal `2x2` untitary of element type `Complex{T}`.
"""
function random_diagonal_unitary(::Type{T}=Float64) where T
    diagm(Vector2(cispi(rand(T)), cispi(rand(T))))
end

function LinearAlgebra.isdiag(m::AbstractMatrix2x2)
    iszero(m[2]) && iszero(m[3])
end

function isantidiag(m::AbstractMatrix2x2)
    iszero(m[1]) && iszero(m[4])
end

isSU2(m::AbstractMatrix2x2) = isone(det(m))

function isSU2(m::AbstractMatrix, approx::AbstractApprox)
    isone(det(m), approx)
end

###
### SU2
###

# It would be easier, if redundant to keep u and t

LinearAlgebra.isdiag(U::AbstractSU2) = iszero(unitary_t(U))

struct SU2{T} <: AbstractSU2{T}
    u::T
    t::T
end

struct SU2B{T, V} <: AbstractSU2{Complex{T}}
    uabs2::T
    alpha_u::V
    alpha_t::V
end

function Base.eltype(::Type{<:SU2B{T}}) where T
    Complex{T}
end

function Base.adjoint(U::SU2B)
    SU2B(U.uabs2, -U.alpha_u, U.alpha_t + pi)
end

# gamma is the angle giving the magnitude.
# It is usually called theta.
# But that's also used for z-rotation angle
struct SU2C{T, V, W} <: AbstractSU2{T}
    function SU2C(gamma, alpha_u::W, alpha_t::W) where {W}
        V = typeof(gamma)
        T = float(V)
        new{T, V, W}(gamma, alpha_u, alpha_t)
    end

    gamma::V
    alpha_u::W
    alpha_t::W
end

@inline SU2C(U::SU2B) = SU2C(acos(sqrt(U.uabs2)), U.alpha_u, U.alpha_t)
@inline unitary_u(su::SU2C) = cos(su.gamma) * cis(su.alpha_u)
@inline unitary_t(su::SU2C) = sin(su.gamma) * cis(su.alpha_t)
@inline unitary_abs2t(U::SU2B) = 1 - U.uabs2
@inline unitary_abs2t(U::AbstractUnitary2x2) = abs2(unitary_t(U))
@inline unitary_u(s::SU2B) = sqrt(s.uabs2) * cis(s.alpha_u)
@inline unitary_t(s::SU2B) = sqrt(1 - s.uabs2) * cis(s.alpha_t)
@inline unitary_abs2u(U::SU2C) = cos(U.gamma)^2
@inline unitary_abs2t(U::SU2C) = sin(U.gamma)^2

@inline function SU2B_from_u_t(u, t)
    abs2u = abs2(u)
    abs2u = abs2u < 1 ? abs2u : one(abs2u)
    alpha_u = iszero(abs2u) ? zero(abs2u) : angle(u / sqrt(abs2u))
    alpha_t = isone(abs2u) ? zero(abs2u) : angle(t / sqrt(1-abs2u))
    SU2B(abs2u, alpha_u, alpha_t)
end

@inline Base.:-(a::SU2B, b::AbstractMatrix2x2) = Matrix2x2(a) - b
@inline Base.:-(a::AbstractMatrix2x2, b::SU2B) = b - a
@inline Base.:+(a::SU2B, b::AbstractMatrix2x2) = Matrix2x2(a) + b
@inline Base.:+(a::AbstractMatrix2x2, b::SU2B) = b + a

@inline Base.:-(a::SU2B, b::SU2B) = Matrix2x2(a) - Matrix2x2(b)
@inline Base.:+(a::SU2B, b::SU2B) = Matrix2x2(a) + Matrix2x2(b)

@inline function Base.:*(a::SU2B, b::SU2B)
    # a2 = SU2(a)
    # b2 = SU2(b)
    # SU2(a2 * b2)

    # Following is slightly faster. I don't know why
    ua = unitary_u(a)
    ta = unitary_t(a)
    ub = unitary_u(b)
    tb = unitary_t(b)
    unew = ua * ub - ta' * tb
    newt = ub * ta + tb * ua'
    SU2B_from_u_t(unew, newt)
end

alt_random_unitary2x2() = alt_random_unitary2x2(Float64)

function alt_random_unitary2x2(::Type{T}) where {T <: AbstractFloat}
    uabs2 = rand(T) # cos^2(gamma)
    uabs = sqrt(uabs2)
    tabs = sqrt(1 - uabs2)
    alpha_u = T(2) * rand(T)
    alpha_t = T(2) * rand(T)
    gamma = T(2) * rand(T)
    pu = cispi(alpha_u + gamma)
    pt = cispi(alpha_t + gamma)
    pua = cispi(-alpha_u + gamma)
    pta = cispi(-alpha_t + gamma)
    Matrix2x2(pu * uabs, pt * tabs, -pta * tabs,  pua * uabs)
end

@inline function Matrix2x2(U::AbstractSU2)
    Matrix2x2(SU2(unitary_u(U), unitary_t(U)))
end

struct ZRot{T} <: AbstractSU2{T}
    function ZRot(t::T) where T
        new{T}(t)
    end
    minushalftheta::T
end

function _float_type(::Type{ZRot{T}}) where T
    float(T)
end

function Base.eltype(::Type{V}) where {V<:ZRot}
    Complex{_float_type(V)}
end

# TODO: organize an interface
# function mulphase(phase, zr::ZRot)
#     u = phase * unitary_u(zr)
#     SU2(u, zero(u))
# end

function Base.:*(x::Number, z::ZRot)
    u = unitary_u(z)
    Matrix2x2(x*u, zero(0), zero(0), x * conj(u))
end

Base.:*(z::ZRot, x::Number) = x * z

# This always gives "positive" rotation.
# Tr(z) >= 0
function random_ZRot(::Type{T}=Float64) where {T<:AbstractFloat}
    a = random_angle(T) / 2
    ZRot(a)
end

#function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2})
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ZRot})
    random_ZRot()
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{ZRot{T}}) where T
    random_ZRot(T)
end

@inline tr(zr::ZRot) = 2 * cos(zr.minushalftheta)

unitary_u(z::ZRot) = cis(z.minushalftheta)
unitary_t(z::ZRot{V}) where {V} = zero(Complex{_float_type(ZRot{V})})
#unitary_t(z::ZRot{<:Any, V}) where {V} = zero(Complex{V})

function eigvals(z::ZRot)
    u = unitary_u(z)
    (u, u')
end

# Probably no faster than fallback
LinearAlgebra.isdiag(::ZRot) = true

Base.:(==)(zr1::ZRot, zr2::ZRot) = zr1.minushalftheta == zr2.minushalftheta

function zrot(theta)
    ZRot(-intdiv(theta, 2))
end

function get_theta(rz::ZRot)
   - 2 * rz.minushalftheta
end

function Base.:*(rz1::ZRot, rz2::ZRot)
    ZRot(rz1.minushalftheta + rz2.minushalftheta)
end

function Base.:/(rz1::ZRot, rz2::ZRot)
    ZRot(rz1.minushalftheta - rz2.minushalftheta)
end

function Base.inv(rz::ZRot)
    ZRot(-rz.minushalftheta)
end

function _add_ZRot(a::ZRot, b::ZRot, op)
    ua = unitary_u(a)
    ub = unitary_u(b)
    u = op(ua, ub)
    Matrix2x2(u, zero(u), zero(u), u')
end

Base.:+(a::ZRot, b::ZRot) = _add_ZRot(a, b, +)
Base.:-(a::ZRot, b::ZRot) = _add_ZRot(a, b, -)

Base.isone(rz::ZRot) = iszero(rz.minushalftheta)
Base.one(rz::ZRot) = ZRot(zero(rz.minushalftheta))

Base.adjoint(rz::ZRot) = ZRot(-rz.minushalftheta)

function SU2B(rz::ZRot{T}) where {T}
    SU2B(one(float(T)), rz.minushalftheta, zero(T))
end

function SU2B(m::Matrix2x2)
    (u, t, c, d) = elements(m)
    uabs2 = abs2(u)
    tabs2 = 1 - uabs2
    if !iszero(uabs2)
        alpha_u = angle(u/sqrt(uabs2))
    else
        alpha_u = zero(uabs2)
    end
    if !iszero(tabs2)
        alpha_t = angle(t/sqrt(tabs2))
    else
        alpha_t = zero(tabs2)
    end
    SU2B(uabs2, radtodar(alpha_u), radtodar(alpha_t))
#    SU2B(uabs2, alpha_u, alpha_t)
end

SU2C(m::Matrix2x2) = SU2C(SU2B(m))

Matrix2x2(rz::ZRot) = Matrix2x2(SU2B(rz))

det(m::AbstractSU2{T}) where {T} = one(T)

# This should not be faster than generic. but it is
function tr(u::SU2B)
    (; uabs2, alpha_u, alpha_t) = u
    uabs = sqrt(uabs2)
    2 * uabs * cos(alpha_u)
end

function eigvals(u::SU2B)
    (; uabs2, alpha_u, alpha_t) = u
    uabs = sqrt(uabs2)
    (s, c) = sincos(alpha_u)
    re_u = c * uabs
    im_u2 = s^2 * uabs2
    tabs2 = 1 - uabs2
    discr = sqrt(im_u2 + tabs2)
    (Complex(re_u, discr), Complex(re_u, -discr))
end



Matrix2x2(U::SU2) = Matrix2x2(U.u, U.t, -U.t', U.u')

SU2(u::SU2) = u
SU2(m::AbstractMatrix) = SU2(m[1,1], m[2,1])

@inline unitary_u(U::SU2) = U.u
@inline unitary_t(U::SU2) = U.t
#@inline SU2_from_u_t(u, t) = SU2(u, t)
@inline SU2_alpha_u(U::SU2) = angle(U.u/abs(U.u))
@inline SU2_alpha_t(U::SU2) = angle(U.t/abs(U.t))

@inline function SU2B(U::SU2)
    SU2B(abs2(U.u), SU2_alpha_u(U), SU2_alpha_t(U))
end

@inline function SU2B(U::SU2C)
    SU2B(unitary_abs2u(U), U.alpha_u, U.alpha_t)
end

@inline function SU2(U::AbstractSU2)
    SU2(unitary_u(U), unitary_t(U))
end

@inline function Base.:-(x::SU2, y::SU2)
    u = x.u - y.u
    t = x.t - y.t
    Matrix2x2(u, t, -t', u')
end

@inline function Base.:+(x::SU2, y::SU2)
    u = x.u + y.u
    t = x.t + y.t
    Matrix2x2(u, t, -t', u')
end

@inline function Base.:*(a::SU2, b::SU2)
    SU2(a.u * b.u - a.t' * b.t, a.t * b.u + a.u' * b.t)
end

function SU2(zr::ZRot)
    SU2(unitary_u(zr), unitary_t(zr))
end

Base.adjoint(u::SU2) = SU2(conj(u.u), -u.t)

"""
    Unitary2x2{T, SUT <: AbstractSU2, V} <: AbstractUnitary2x2{T}

Represents a `2 x 2` unitary matrix as a matrix in SU2 times a global phase.

Fields: `su2`, `phi`.
"""
struct Unitary2x2{T, SUT <: AbstractSU2, V} <: AbstractUnitary2x2{T}
    function Unitary2x2(su2::W, phi::V) where {V, W <: AbstractSU2{T}} where {T}
        new{T, W, V}(su2, phi)
    end

    su2::SUT
    phi::V
end

function Base.adjoint(u::Unitary2x2)
    Unitary2x2(adjoint(u.su2), -u.phi)
end

function Unitary2x2(su::AbstractSU2{T}) where T
    Unitary2x2(su, zero(T))
end

function elements(m::Unitary2x2)
    ph = cis(m.phi)
    map(x -> ph * x, elements(m.su2))
end

det(U::Unitary2x2) = cis(2 * U.phi)

function Unitary2x2(m::Matrix2x2{Complex{T}}, ::Type{SU2T} = SU2) where {T <: AbstractFloat, SU2T <: AbstractSU2}
    phase_fac = sqrt(det(m))
    msu2 = map(x -> x / phase_fac, m)
    Unitary2x2(SU2T(msu2), radtodar(angle(phase_fac)))
end

function eigvals(U::Unitary2x2)
    (; su2, phi) = U
    (v1, v2) = eigvals(su2)
    p = cis(phi)
    (p * v1, p * v2)
end

function Matrix2x2(U::Unitary2x2)
    su2 = U.su2
    phi = U.phi
    p = cis(phi)
    um = Matrix2x2(su2)
    p * um
end

random_unitary2x2() = random_unitary2x2(Float64)

function random_unitary2x2(::Type{T}) where {T <: AbstractFloat}
    su2 = random_SU2B(T)
    gamma = T(2) * rand(T) # Huh random_angle
    Unitary2x2(su2, Dar(gamma))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Matrix2x2})
    rand(rng, Matrix2x2{Float64})
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Matrix4x4})
    rand(rng, Matrix4x4{Float64})
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{MatrixNxN{T, N, M}}) where {T, N, M}
    # Using ntuple is very slow and allocates. Even with compile time length
    # tup = ntuple(_  -> rand(T), M)
    # This is a bit faster, even though it allocates
    tup = NTuple{M, T}(rand(T, M))
    MatrixNxN{T, N, M}(tup)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{MatrixNxN{T, 4, 16}}) where {T}
    # ntuple is very slow. allocating with rand(16) is less slow, but still slow
    #    tup = ntuple(_  -> rand(T), 16)
    #    MatrixNxN(NTuple{16, T}(tup))
    # Writing rand(T) sixteen times is efficient. Maybe use st like a generated function?
    tup = (rand(T),rand(T),rand(T),rand(T),
           rand(T),rand(T),rand(T),rand(T),
           rand(T),rand(T),rand(T),rand(T),
           rand(T),rand(T),rand(T),rand(T))
    MatrixNxN{T, 4, 16}(tup)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{MatrixNxN{T, 2, 4}}) where {T}
    # Both ntuple and unrolled are equally performant
    # tup = ntuple(_  -> rand(T), 4)
    tup = (rand(T),rand(T),rand(T),rand(T))
    MatrixNxN{T, 2, 4}(NTuple{4, T}(tup))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2})
    rand(rng, Random.SamplerType{SU2{ComplexF64}}())
end

# Broken
function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2{Complex{T}}}) where {T<:Real}
#    samp = Random.SamplerType{SU2B{Complex{T}}}()
    samp = Random.SamplerType{SU2B{T}}()
    inn = rand(rng, samp)
    SU2(inn)
#    SU2(rand(rng, Random.SamplerType{SU2B{Complex{T}}}()))
end

function random_SU2B(rng::Random.AbstractRNG)
    random_SU2B(rng, Float64)
end

function random_SU2B()
    random_SU2B(Float64)
end

function random_SU2B(::Type{FloatT}) where {FloatT}
    random_SU2B(Random.GLOBAL_RNG, FloatT)
end

# Julia's random sampling interface is a bit difficult. We do this manually.
function random_SU2B(rng::Random.AbstractRNG, ::Type{FloatT}) where {FloatT}
    T = Float64
    uabs2 = rand(rng, T) # cos^2(gamma)
    SU2B(uabs2, rand(rng, Dar{T}), rand(rng, Dar{T}))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2B})
    T = Float64
    uabs2 = rand(rng, T) # cos^2(gamma)
    SU2B(uabs2, rand(rng, Dar{T}), rand(rng, Dar{T}))
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2B{T, Dar{T}}}) where {T}
    uabs2 = rand(rng, T) # cos^2(gamma)
    SU2B(uabs2, rand(rng, Dar{T}), rand(rng, Dar{T}))
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2B{Complex{T}}}) where {T}
    uabs2 = rand(rng, T) # cos^2(gamma)
    SU2B(uabs2, rand(rng, Ang), rand(rng, Ang))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2C})
    rand(rng, Random.SamplerType{SU2C{ComplexF64}}())
end

# BROKEN
function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2C{Complex{T}}}) where {T}
    SU2C(rand(rng, Random.SamplerType{SU2B{Complex{T}}}()))
end

# These all return `Ang`. But that only handles fixed precision
# what about BigFloat
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Unitary2x2{T}}) where T
    Unitary2x2(rand(rng, SU2B{Complex{T}}), rand(rng, Ang))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Unitary2x2})
    Unitary2x2(rand(rng, SU2B), rand(rng, Ang))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Unitary2x2{T, SUT, V}}) where {V, SUT <: AbstractSU2{T}} where {T}
    Unitary2x2(rand(rng, SUT), rand(rng, Ang))
end

struct ScaleMatrix2x2{V, MatrixT <: AbstractMatrix2x2, ScaleT} <: AbstractMatrix2x2{V}
    function ScaleMatrix2x2(m::AbstractMatrix2x2, s)
        if isa(m, ScaleMatrix2x2)
            throw(ArgumentError(lazy"Nesting ScaleMatrix not supported."))
        end
        V = typeof(s * m[1])
        new{V, typeof(m), typeof(s)}(m, s)
    end

    m::MatrixT
    s::ScaleT
end

function scalematrix end

for f in (:real, :imag, :float, :big, :complex, :canonical)
    @eval ($f)(m::Matrix2x2) = map($f, m)
end

for f in (:real, :imag, :float, :big, :complex)
    @eval ($f)(m::ScaleMatrix2x2) = map($f, Matrix2x2(m))
end


function Base.getindex(sm::ScaleMatrix2x2, i::Integer)
    sm.s * sm.m[i]
end

function Matrix2x2(sm::ScaleMatrix2x2)
    Matrix2x2(map(x -> sm.s * x, sm.m))
end

function Base.map(f, sm::ScaleMatrix2x2)
    ScaleMatrix2x2(map(f, sm.m), sm.s)
end

function Base.show(io::IO, ::PRETTY, sm::ScaleMatrix2x2)
    summary(io, sm)
    println(io, ":")
    show(io, PRETTY(), sm.s)
    println(io, " ×")
    _show_matrix2x2(io, sm.m)
end

function LinearAlgebra.eigvals(sm::ScaleMatrix2x2)
    (v1, v2) = eigvals(sm.m)
    (sm.s * v1, sm.s * v2)
end

function Base.adjoint(sm::ScaleMatrix2x2)
    ScaleMatrix2x2(adjoint(sm.m), adjoint(sm.s))
end

function Base.:*(sm1::ScaleMatrix2x2, sm2::ScaleMatrix2x2)
    s = sm1.s *  sm2.s
    m = sm1.m * sm2.m
    ScaleMatrix2x2(m, s)
end

function Base.:*(m1::Matrix2x2, sm2::ScaleMatrix2x2)
    sm1 = scalematrix(m1)
    s = sm1.s *  sm2.s
    m = sm1.m * sm2.m
    ScaleMatrix2x2(m, s)
end
Base.:*(sm2::ScaleMatrix2x2, m1::Matrix2x2) = m1 * sm2

# Fallback in dense.jl is no less performant in at least some cases.
function Base.:^(sm::ScaleMatrix2x2, n::Integer)
    ScaleMatrix2x2(sm.m^n, sm.s^n)
end

det(sm::ScaleMatrix2x2) = (sm.s * sm.s) * det(sm.m)

end # module Matrices2x2
