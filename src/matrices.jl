module Matrices2x2

import Random
import LinearAlgebra: eigvals, svdvals, opnorm, tr, det, diag, diagm, eigvecs, norm, normalize,
    dot
import LinearAlgebra
import ..Common: canonical
import ..Utils: PRETTY, cpad, _show_with_fieldnames, _power_by_squaring
import ..Angles: radtodar, Dar, scalepi, unscalepi, intdiv
import IsApprox: isunitary, AbstractApprox, Equal, Approx

abstract type AbstractMatrix2x2{T} <: AbstractMatrix{T} end
abstract type Normal2x2{T} <: AbstractMatrix2x2{T} end
abstract type AbstractUnitary2x2{T} <: Normal2x2{T} end
abstract type AbstractSU2{T} <: AbstractUnitary2x2{T} end

abstract type AbstractVector2{T} <: AbstractVector{T} end

"""
    Matrix2x2{T} <: AbstractMatrix{T}

Stack allocated, immutable, 2 x 2 matrix.

If `T` is `isbits`, then `Matrix2x2{T}` is `isbits`.

We implement a few necessary operations, including matrix multiplication,
addition, subtraction, and unary minus.
"""
struct Matrix2x2{T} <: AbstractMatrix2x2{T}
    data::NTuple{4, T}
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
    size(m) == (2,2) || error("Wrong size")
    Matrix2x2(m...,)
end

Base.Matrix(m::AbstractMatrix2x2) = reshape([elements(m)...,], (2,2))

Matrix2x2{T}(m::AbstractMatrix2x2) where {T} = map(x -> T(x), elements(m))

Matrix2x2(m::Matrix2x2) = m

function Base.map(f, m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    Matrix2x2(f(a), f(b), f(c), f(d))
end

elements(m::Matrix2x2) = m.data
elements(m::AbstractMatrix2x2) = elements(Matrix2x2(m))

##
## Conversion and related construction
##

Base.convert(::Type{Matrix2x2{T}}, m::Matrix2x2) where {T} = Matrix2x2{T}(m)
Base.AbstractFloat(m::Matrix2x2) = map(float, m)
#Base.float(m::Matrix2x2) = AbstractFloat(m)
#Base.complex(m::Matrix2x2) = map(complex, m)
# This is not conventional. Change this
# Base.big(m::Matrix2x2) = map(big, m)
# Base.real(m::Matrix2x2) = map(real, m)
# Base.imag(m::Matrix2x2) = map(imag, m)

import Base: real, imag, big, complex, float

for f in (:real, :imag, :float, :big, :complex, :canonical)
    @eval ($f)(m::Matrix2x2) = map($f, m)
end


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
    spc = "  "
    (as, bs, cs, ds)  = map(_showstr, elements(m))
    (al, bl, cl, dl) = map(length, (as, bs, cs, ds))
    w1 = max(al, bl)
    w2 = max(cl, dl)
    print(io, cpad(as, w1), spc)
    println(io, cpad(cs, w2))
    print(io, cpad(bs, w1), spc)
    print(io, cpad(ds, w2))
end

function Base.show(io::IO, m::Matrix2x2)
    print(io, typeof(m), "(")
    i = 0
    for x in element(m)
        i += 1
        if i > 1
            print(io, ", ")
        end
        print(io, x)
    end
    print(io, ")")
end

function Base.show(io::IO, ::PRETTY, v::AbstractVector2)
    (v1, v2) = elements(v)
#    summary(io, v)
#    println(io, ":")
    show(io, PRETTY(), v1)
#    println(io)
    println(io)
    show(io, PRETTY(), v2)
end


# function Base.show(io::IO, ::PRETTY, m::Matrix2x2)
#     summary(io, m)
#     println(io, ":")
#     spc = "  "
#     (as, bs, cs, ds)  = map(_showstr, m.data)
#     (al, bl, cl, dl) = map(length, (as, bs, cs, ds))
#     w1 = max(al, bl)
#     w2 = max(cl, dl)
#     print(io, cpad(as, w1), spc)
#     println(io, cpad(cs, w2))
#     print(io, cpad(bs, w1), spc)
#     print(io, cpad(ds, w2))
# end

##
## Required, and standard, `Base` properties
##

@inline function Base.getindex(m::Matrix2x2, i::Integer)
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

Base.size(::AbstractMatrix2x2) = (2, 2)
Base.eltype(::Type{<:AbstractMatrix2x2{T}}) where T = T
Base.eltype(::Type{AbstractMatrix2x2{T}}) where T = T
Base.IndexStyle(::Type{<:AbstractMatrix2x2}) = IndexLinear()
Base.one(::Type{Matrix2x2{T}}) where T = Matrix2x2(one(T), zero(T), zero(T), one(T))
Base.one(::Matrix2x2{T}) where {T} = one(Matrix2x2{T})
Base.zero(::Type{Matrix2x2{T}}) where T = Matrix2x2(zero(T), zero(T), zero(T), zero(T))
Base.zero(::Matrix2x2{T}) where {T} = zero(Matrix2x2{T})

function Base.isone(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    isone(a) && isone(d) && iszero(b) && iszero(c)
end

function Base.iszero(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    iszero(a) && iszero(d) && iszero(b) && iszero(c)
end

##
## Further properties.
##

##
## IsApprox.isunitary
##

# function isunitary(m::Matrix2x2)
#     (a, b, c, d) = m.data
#     (abs2(a) + abs2(c) == one(a)) &&
#         (abs2(b) + abs2(d) == one(a)) &&
#         ((conj(a) * b + d * conj(c)) == zero(a))
# end

function isunitary(m::Matrix2x2)
    isone(m * m')
end


isunitary(m::Matrix2x2, ::Equal) = isunitary(m)

isunitary(m::AbstractUnitary2x2) = true
svdvals(::AbstractUnitary2x2{T}) where {T} = (one(real(T)), one(real(T)))

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

@inline function Base.:*(m::Matrix2x2, x::Number)
    (a, b, c, d) = elements(m)
    Matrix2x2(x * a, x * b, x * c, x * d)
end

@inline Base.:*(x::Number, m::Matrix2x2) = m * x

@inline function _pair_op(op, m1, m2)
    (a1, b1, c1, d1) = elements(m1)
    (a2, b2, c2, d2) = elements(m2)
    Matrix2x2(op(a1, a2), op(b1, b2), op(c1, c2), op(d1,d2))
end

Base.:*(m1::Matrix2x2, m2::Matrix2x2) = _mul2(m1, m2)
Base.:+(m1::Matrix2x2, m2::Matrix2x2) = _pair_op(+, m1, m2)
Base.:-(m1::Matrix2x2, m2::Matrix2x2) = _pair_op(-, m1, m2)
Base.:-(m::Matrix2x2) = map(-, m) # Unary minus

Base.:-(a::AbstractMatrix2x2, b::AbstractMatrix2x2) = Matrix2x2(a) - Matrix2x2(b)

# This is only called if `n` is not literal at the call site.
function Base.:^(m::Matrix2x2, n::Integer)
    n == 0 && return one(m)
    n == 1 && return m  # If elements of m are mutable, this is problematic
    n == 2 && return m * m
    return _power_by_squaring(m, n)
end

Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{0}) = one(m)
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{1}) = m # If elements of m are mutable, this is problematic
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{2}) = m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{3}) = m * m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{4}) = (m * m) * (m * m)

##
## Linear algebra and related operations
##

function Base.permutedims(m::Matrix2x2)
    (a, b, c, d) = elements(m)
    Matrix2x2(a, c, b, d)
end

Base.adjoint(m::Matrix2x2) = permutedims(map(adjoint, m))
Base.transpose(m::Matrix2x2) = permutedims(map(transpose, m))

@inline function tr(m::AbstractMatrix2x2)
    (a, _b, _c, d) = elements(m)
    a + d
end

@inline function det(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    canonical(a * d - c * b)
end

@inline tr(m::AbstractSU2) = 2*real(m[1])

"""
    eigvals(m::Matrix2x2)::NTuple{2}

Return a tuple of the eigenvalues of `m`.
"""
function eigvals(m::AbstractMatrix2x2)
    (a, b, c, d) = elements(m)
    discr = sqrt((a-d)^2 + 4*b*c)
    (
        (a + d + discr)/2,
        (a + d - discr)/2,
    )
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

function eigvecs(m::AbstractMatrix2x2)
    (v1, v2) = eigvals(m)
    (a, b, c, d) = elements(m)
    vec1 = normalize(Vector2(c, v1 - a))
    vec2 = normalize(Vector2(v2 - d, b))
    Matrix2x2(vec1[1], vec1[2], vec2[1], vec2[2])
end

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
function svdvals(m::Matrix2x2)
    ma = m * adjoint(m)
    if isone(ma)
        ev = one(real(eltype(m)))
        return (ev, ev)
    end
    (v1, v2) = eigvals(ma)
    (s1, s2) = (sqrt(real(v1)), sqrt(real(v2)))
    s1 > s2 ? (s1,s2) : (s2, s1)
end

function opnorm(m::AbstractMatrix2x2)
    (v1, v2) = svdvals(m)
    max(v1, v2)
end

function opnormdistance(a::AbstractMatrix2x2, b::AbstractMatrix2x2)
    opnorm(a - b)
end

"""
    tracenorm(m::Matrix2x2)

Return the trace norm of `m`.
"""
function tracenorm(m::AbstractMatrix2x2)
    (v1, v2) = svdvals(m)
    v1 + v2
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

###
### Parameterizations of unitaries
###

struct UnitaryParam1{T}
    u::Complex{T}
    t::Complex{T}
    phi::T
end

struct UnitaryParam2{T}
    gamma::T
    alpha_u::T
    alpha_t::T
    phi::T
end

function Base.show(io::IO, ::PRETTY, U::UnitaryParam1)
    _show_with_fieldnames(io, U)
end

function Base.show(io::IO, ::PRETTY, U::UnitaryParam2)
    _show_with_fieldnames(io, U)
end

function unitary_compose(U::UnitaryParam1)
    (;u, t, phi) = U
    cp = cis(phi)
    Matrix2x2(u, t, -conj(t)*cp, conj(u)*cp)
end

function unitary_compose(U::UnitaryParam2)
    (;gamma, alpha_u, alpha_t, phi) = U
    u = cis(alpha_u) * cos(gamma)
    t = cis(alpha_t) * sin(gamma)
    unitary_compose(UnitaryParam1(u, t, phi))
end

function unitary_decompose(U::Matrix2x2, ::Type{T}) where {T<:UnitaryParam1}
    (u, t, c, d) = elements(U)
    if !iszero(u)
        ua = abs2(u)
        u_d = u * d # |u|^2 cis(phi)
        cp = u_d / ua # |u|^2/|u|^2 cis(phi) == cis(phi)
        phi = angle(cp) # angle(cis(phi)) == phi
    else
        phi = zero(u)
    end
    return T(u, t, phi)
end

function unitary_decompose(U::Matrix2x2, ::Type{T}) where {T<:UnitaryParam2}
    (;u, t, phi) = unitary_decompose(U, UnitaryParam1)
    (gamma, alpha_u, alpha_t) = _unitary_decompose_special(u, t)
    return T(gamma, alpha_u, alpha_t, phi)
end

function _unitary_decompose_special(u, t)
    uabs = abs(u)
    tabs = abs(t)
    if uabs > tabs
        gamma = acos(uabs)
    else
        gamma = asin(tabs)
    end
    alpha_u = iszero(uabs) ? zero(uabs) : angle(u/uabs)
    alpha_t = iszero(tabs) ? zero(tabs) : angle(t/tabs)
    return (gamma, alpha_u, alpha_t)
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

struct SU2B{T, V} <: AbstractSU2{T}
    uabs2::T
    alpha_u::V
    alpha_t::V
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

@inline function SU2B_from_u_t(u, t)
    abs2u = abs2(u)
    SU2B(abs2u, angle(u/sqrt(abs2u)), angle(t/sqrt(1-abs2u)))
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

random_unitary2x2() = random_unitary2x2(Float64)

function random_unitary2x2(::Type{T}) where {T <: AbstractFloat}
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

struct ZRot{T, V} <: AbstractSU2{T}
    function ZRot(t::T) where T
        V = float(T)
        new{T, V}(t)
    end
    minushalftheta::T
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
function random_ZRot()
    a = random_angle() / 2
    ZRot(a)
end

@inline tr(zr::ZRot) = 2 * cos(zr.minushalftheta)

unitary_u(z::ZRot) = cis(z.minushalftheta)
unitary_t(z::ZRot{<:Any, V}) where {V} = zero(Complex{V})

function eigvals(z::ZRot)
    u = unitary_u(z)
    (u, u')
end

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

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2})
    rand(rng, Random.SamplerType{SU2{ComplexF64}}())
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2{Complex{T}}}) where {T}
    SU2(rand(rng, Random.SamplerType{SU2B{Complex{T}}}()))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2B})
    rand(rng, Random.SamplerType{SU2B{ComplexF64}}())
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2B{Complex{T}}}) where {T}
    uabs2 = rand(T) # cos^2(gamma)
    alpha_u = T(2) * rand(T)
    alpha_t = T(2) * rand(T)
    SU2B(uabs2, Dar(alpha_u), Dar(alpha_t))
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{SU2C})
    rand(rng, Random.SamplerType{SU2C{ComplexF64}}())
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerType{SU2C{Complex{T}}}) where {T}
    SU2C(rand(rng, Random.SamplerType{SU2B{Complex{T}}}()))
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

struct Unitary2x2{T, SUT <: AbstractSU2, V} <: AbstractUnitary2x2{T}
    function Unitary2x2(su2::W, phi::V) where {V, W <: AbstractSU2{T}} where {T}
        new{T, W, V}(su2, phi)
    end

    su2::SUT
    phi::V
end

det(U::Unitary2x2) = cis(2 * U.phi)

function Unitary2x2(su::AbstractSU2{T}) where T
    Unitary2x2(su, zero(T))
end

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
    (;su2, phi) = U
    p = cis(phi)
    map(x -> p * x, Matrix2x2(su2))
end

end # module Matrices2x2
