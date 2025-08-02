module Matrices
import LinearAlgebra: eigvals, svdvals, opnorm, tr, det
import ..Common: canonical
import ..Utils: PRETTY, cpad

"""
    Matrix2x2{T} <: AbstractMatrix{T}

Stack allocated, immutable 2 x 2 matrix.

If `T` is `isbits`, then `Matrix2x2{T}` is `isbits`.

We implement a few necessary operations: matrix multiplication,
addition, subtraction, unary minus. Most other things will be computed via
fallback methods, which may be inefficient. For example eigenvalues. We could implement
these if necessary.
"""
struct Matrix2x2{T} <: AbstractMatrix{T}
    data::NTuple{4, T}
end

# Hmmm. Depends on whether type T is copyable
# In any case, a Tuple is not.
# Base.copy(m::Matrix2x2) = Matrix2x2(copy(m.data))

Matrix2x2(a, b, c, d) = Matrix2x2(promote(a, b, c, d))
Matrix2x2{T}(a, b, c, d) where {T} = Matrix2x2(T(a), T(b), T(c), T(d))

Base.size(::Matrix2x2) = (2, 2)
Base.eltype(::Type{Matrix2x2{T}}) where T = T
Base.IndexStyle(::Type{<:Matrix2x2}) = IndexLinear()
Base.one(::Type{Matrix2x2{T}}) where T = Matrix2x2(one(T), zero(T), zero(T), one(T))
Base.one(::Matrix2x2{T}) where {T} = one(Matrix2x2{T})
Base.zero(::Type{Matrix2x2{T}}) where T = Matrix2x2(zero(T), zero(T), zero(T), zero(T))
Base.zero(::Matrix2x2{T}) where {T} = zero(Matrix2x2{T})

function _showstr(obj)
    b = IOBuffer()
    show(IOContext(b, :compact=>true), PRETTY(), obj)
    return String(take!(b))
end

function Base.show(io::IO, ::PRETTY, m::Matrix2x2)
    summary(io, m)
    println(io, ":")
    spc = "  "
    (as, bs, cs, ds)  = map(_showstr, m.data)
    (al, bl, cl, dl) = map(length, (as, bs, cs, ds))
    w1 = max(al, bl)
    w2 = max(cl, dl)
    print(io, cpad(as, w1), spc)
    println(io, cpad(cs, w2))
    print(io, cpad(bs, w1), spc)
    print(io, cpad(ds, w2))
end

function Base.map(f, m::Matrix2x2)
    (a, b, c, d) = m.data
    Matrix2x2(f(a), f(b), f(c), f(d))
end

@inline function Base.getindex(m::Matrix2x2, i::Integer)
    @boundscheck checkbounds(m, i)
    return @inbounds m.data[i]
end

@inline function _mul2(m1, m2)
    (a1, b1, c1, d1) = m1.data
    (a2, b2, c2, d2) = m2.data
    Matrix2x2(a1*a2 + c1*b2, a2*b1 + d1*b2, a1*c2 + c1*d2, b1*c2 + d1*d2)
#    Matrix2x2(a1*a2 + c1*b2, a1*c2+c1*d2, a2*b1+d1*b2, b1*c2 + d1*d2)
#    Matrix2x2(a1*a2 + b1*c2, a1*b2+b1*d2, a2*c1+d1*c2, c1*b2 + d1*d2)
end

@inline function Base.:*(m::Matrix2x2, x::Number)
    (a, b, c, d) = m.data
    Matrix2x2(x * a, x * b, x * c, x * d)
end

@inline Base.:*(x::Number, m::Matrix2x2) = m * x

@inline function _pair_op(op, m1, m2)
    (a1, b1, c1, d1) = m1.data
    (a2, b2, c2, d2) = m2.data
    Matrix2x2(op(a1, a2), op(b1, b2), op(c1, c2), op(d1,d2))
end

Base.:*(m1::Matrix2x2, m2::Matrix2x2) = _mul2(m1, m2)
Base.:+(m1::Matrix2x2, m2::Matrix2x2) = _pair_op(+, m1, m2)
Base.:-(m1::Matrix2x2, m2::Matrix2x2) = _pair_op(-, m1, m2)
Base.:-(m::Matrix2x2) = map(-, m) # Unary minus

# This is only called if `n` is not literal at the call site.
function Base.:^(m::Matrix2x2, n::Integer)
    n == 0 && return one(m)
    n == 1 && return m  # If elements of m are mutable, this is problematic
    n == 2 && return m * m
    return Base.power_by_squaring(m, n)
end

Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{0}) = one(m)
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{1}) = m # If elements of m are mutable, this is problematic
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{2}) = m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{3}) = m * m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{4}) = (m * m) * (m * m)

Matrix2x2{T}(m::Matrix2x2) where {T} = map(x -> convert(T, x), m)
Base.convert(::Type{Matrix2x2{T}}, m::Matrix2x2) where {T} = Matrix2x2{T}(m)
Base.float(m::Matrix2x2) = AbstractFloat(m)
Base.complex(m::Matrix2x2) = float(m)
# This is not conventional. Change this
Base.AbstractFloat(m::Matrix2x2) = map(float, m)
Base.big(m::Matrix2x2) = map(big, m)

"""
    canonical(m::Matrix2x2)

Return a new matrix by calling `canonical` element-wise on `m`.
"""
canonical(m::Matrix2x2) = map(canonical, m)

"""
    get_theta(m::Matrix2x2)

Find `theta` from a Z-rotation matrix `m` with possible global phase.

Assume `m` is diagonal, with `m[1] = cis(-theta/2 + phi)`
and `m[4] = cis(theta/2 + phi)`. Return `theta`.
"""
get_theta(m::Matrix2x2) = angle(m[4] / m[1])

function Base.transpose(m::Matrix2x2)
    (a, b, c, d) = m.data
    Matrix2x2(a, c, b, d)
end

function Base.adjoint(m::Matrix2x2)
    (a, b, c, d) = map(adjoint, m.data)
    Matrix2x2(a, c, b, d)
end

tr(m::Matrix2x2) = m[1] + m[4]

det(m::Matrix2x2) = m[1] * m[4] - m[3] * m[2]

function eigvals(m::Matrix2x2)
    (a, b, c, d) = m.data
    discr = sqrt(a^2 + 4*b*c - 2*a*d + d^2)
    (
        (a + d - discr)/2,
        (a + d + discr)/2
    )
end

function svdvals(m::Matrix2x2)
    ma = m * adjoint(m)
    (v1, v2) = eigvals(ma)
    (sqrt(real(v1)), sqrt(real(v2)))
end

function opnorm(m::Matrix2x2)
    (v1, v2) = svdvals(m)
    max(v1, v2)
end


end # module Matrices
