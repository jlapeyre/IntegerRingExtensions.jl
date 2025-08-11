module Matrices2x2
import LinearAlgebra: eigvals, svdvals, opnorm, tr, det, diag, diagm
import ..Common: canonical
import ..Utils: PRETTY, cpad, _show_with_fieldnames
import IsApprox: isunitary, AbstractApprox, Equal, Approx

abstract type AbstractMatrix2x2{T} <: AbstractMatrix{T} end
abstract type Normal2x2{T} <: AbstractMatrix2x2{T} end

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
Base.float(m::Matrix2x2) = AbstractFloat(m)
Base.complex(m::Matrix2x2) = float(m)
# This is not conventional. Change this
Base.AbstractFloat(m::Matrix2x2) = map(float, m)
Base.big(m::Matrix2x2) = map(big, m)

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

Base.size(::AbstractMatrix2x2) = (2, 2)
Base.eltype(::Type{<:AbstractMatrix2x2{T}}) where T = T
Base.eltype(::Type{AbstractMatrix2x2{T}}) where T = T
Base.IndexStyle(::Type{<:AbstractMatrix2x2}) = IndexLinear()
Base.one(::Type{Matrix2x2{T}}) where T = Matrix2x2(one(T), zero(T), zero(T), one(T))
Base.one(::Matrix2x2{T}) where {T} = one(Matrix2x2{T})
Base.zero(::Type{Matrix2x2{T}}) where T = Matrix2x2(zero(T), zero(T), zero(T), zero(T))
Base.zero(::Matrix2x2{T}) where {T} = zero(Matrix2x2{T})

function Base.isone(m::Matrix2x2)
    (a, b, c, d) = m.data
    isone(a) && isone(d) && iszero(b) && iszero(c)
end

function Base.iszero(m::Matrix2x2)
    (a, b, c, d) = m.data
    iszero(a) && iszero(d) && iszero(b) && iszero(c)
end

##
## Further properties.
##

##
## IsApprox.isunitary
##

function isunitary(m::Matrix2x2)
    (a, b, c, d) = m.data
    (abs2(a) + abs2(c) == one(a)) &&
        (abs2(b) + abs2(d) == one(a)) &&
        ((conj(a) * b + d * conj(c)) == zero(a))
end

isunitary(m::Matrix2x2, ::Equal) = isunitary(m)

abstract type AbstractUnitary2x2{T} <: Normal2x2{T} end

isunitary(m::AbstractUnitary2x2) = true
svdvals(::AbstractUnitary2x2{T}) where {T} = (one(T), one(T))

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
    (a, b, c, d) = m.data
    isapprox(abs2(a) + abs2(c), one(a); app.kw...) || return false
    isapprox(abs2(b) + abs2(d), one(a); app.kw...) || return false

    isapprox(conj(a) * b, - d * conj(c); app.kw...)
end

"""
    canonical(m::Matrix2x2)

Return a new matrix by calling `canonical` element-wise on `m`.

This may reduce or canonicalize elements that implement `canonical`.
"""
canonical(m::Matrix2x2) = map(canonical, m)

##
## Matrix arithmetic
##

@inline function _mul2(m1, m2)
    (a1, b1, c1, d1) = m1.data
    (a2, b2, c2, d2) = m2.data
    Matrix2x2(a1*a2 + c1*b2, a2*b1 + d1*b2, a1*c2 + c1*d2, b1*c2 + d1*d2)
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

##
## Linear algebra and related operations
##

function Base.permutedims(m::Matrix2x2)
    (a, b, c, d) = m.data
    Matrix2x2(a, c, b, d)
end

Base.adjoint(m::Matrix2x2) = permutedims(map(adjoint, m))
Base.transpose(m::Matrix2x2) = permutedims(map(transpose, m))
tr(m::Matrix2x2) = m[1] + m[4]
det(m::Matrix2x2) = m[1] * m[4] - m[3] * m[2]

"""
    eigvals(m::Matrix2x2)::NTuple{2}

Return a tuple of the eigenvalues of `m`.
"""
function eigvals(m::Matrix2x2)
    (a, b, c, d) = m.data
    discr = sqrt(a^2 + 4*b*c - 2*a*d + d^2)
    (
        (a + d + discr)/2,
        (a + d - discr)/2,
    )
end

"""
    svdvals(m::Matrix2x2)::NTuple{2}

Return a tuple of the singular values of `m` in descending order.
"""
function svdvals(m::Matrix2x2)
    ma = m * adjoint(m)
    (v1, v2) = eigvals(ma)
    (s1, s2) = (sqrt(real(v1)), sqrt(real(v2)))
    s1 > s2 ? (s1,s2) : (s2, s1)
end

function opnorm(m::Matrix2x2)
    (v1, v2) = svdvals(m)
    max(v1, v2)
end

"""
    tracenorm(m::Matrix2x2)

Return the trace norm of `m`.
"""
function tracenorm(m::Matrix2x2)
    (v1, v2) = svdvals(m)
    v1 + v2
end

"""
    GPID(A::Matrix2x2, B::Matrix2x2)

Compute the global phase invariant distance between `A` and `B`.

(See Mukhopadhyay 2021)
"""
function GPID(A::Matrix2x2, B::Matrix2x2)
    (a, b, c, d) = map(conj, A.data)
    (w, x, y, z) = B.data
    trprod = a*w + b*x + c*y + d*z
    return sqrt(1 - abs(trprod) / 2)
end

##
## Vector2
##

struct Vector2{T} <: AbstractVector{T}
    data::NTuple{2, T}
end

Vector2(a, b) = Vector2(promote(a, b))
Vector2{T}(a, b) where {T} = Vector2(T(a), T(b))
Base.size(::Vector2) = (2,)
Base.eltype(::Type{Vector2{T}}) where T = T
Base.IndexStyle(::Type{<:Vector2}) = IndexLinear()

@inline function Base.getindex(m::Vector2, i::Integer)
    @boundscheck checkbounds(m, i)
    return @inbounds m.data[i]
end

function Base.:*(m::Matrix2x2, v::Vector2)
    (a,b,c,d) = m.data
    (x, y) = v.data
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
    (u, t, c, d) = U.data
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

random_SU2() = random_SU2(Float64)

function random_SU2(::Type{T}) where {T}
    uabs2 = rand(T) # cos^2(gamma)
    alpha_u_scaled = T(2) * rand(T)
    alpha_t_scaled = T(2) * rand(T)
    return SU2(uabs2, alpha_u_scaled, alpha_t_scaled)
end

random_unitary2x2() = random_unitary2x2(ComplexF64)
function random_unitary2x2(::Type{T}) where {T <: AbstractFloat}
    random_unitary2x2(Complex{T})
end

function random_unitary2x2(::Type{Complex{T}}) where {T <: Real}
    su2 = random_SU2(T)
    phi = T(2) * rand(T)
    return Unitary2x2(su2, phi)
    # phi = T(2) * rand(T)
    # (; uabs2, alpha_u, alpha_t) = random_SU2(T)
    # tabs = sqrt(one(uabs2) - uabs2)
    # uabs = sqrt(uabs2)
    # u = cispi(phi + alpha_u) * uabs
    # t = cispi(phi + alpha_t) * tabs
    # c = -cispi(phi - alpha_t) * tabs
    # d = cispi(phi - alpha_u) * uabs
    # return Matrix2x2(u, t, c, d)
end

# function random_SU2(::Type{T}) where {T <: AbstractFloat}
#     random_SU2(Complex{T})
# end

# random_SU2() = random_SU2(ComplexF64)
# function random_SU2(::Type{Complex{T}}) where {T <: Real}
#     (; uabs2, alpha_u_scaled, alpha_t_scaled) = _random_SU2(T)
#     uabs = sqrt(uabs2)
#     tabs = sqrt(one(uabs2) - uabs2)
#     uphase = cispi(alpha_u_scaled)
#     tphase = cispi(alpha_t_scaled)
#     u = uphase * uabs
#     t = tphase * tabs
#     c = -conj(tphase) * tabs
#     d = conj(uphase) * uabs
#     Matrix2x2(u, t, c, d)
# end

struct SU2{T} <: AbstractUnitary2x2{T}
    uabs2::T
    alpha_u::T
    alpha_t::T
end

det(m::SU2{T}) where {T} = one(T)

function eigvals(u::SU2)
    (; uabs2, alpha_u, alpha_t) = u
    uabs = sqrt(uabs2)
    (s, c) = sincospi(alpha_u)
    re_u = c * uabs
    im_u2 = s^2 * uabs2
    tabs2 = 1 - uabs2
    discr = sqrt(im_u2 + tabs2)
    (Complex(re_u, discr), Complex(re_u, -discr))
end

function Matrix2x2(su2::SU2)
    (; uabs2, alpha_u, alpha_t) = su2
    uabs = sqrt(uabs2)
    tabs = sqrt(one(uabs2) - uabs2)
    pu = cispi(alpha_u)
    pt = cispi(alpha_t)
    Matrix2x2(uabs * pu, tabs * pt, - tabs * conj(pt), uabs * conj(pu))
end

struct Unitary2x2{T} <: AbstractUnitary2x2{T}
    su2::SU2{T}
    phi::T
end

function eigvals(U::Unitary2x2)
    (; su2, phi) = U
    (v1, v2) = eigvals(su2)
    p = cispi(phi)
    # Relatively inefficient to do ths.
    (p * v1, p * v2)
end

function Matrix2x2(U::Unitary2x2)
    (;su2, phi) = U
    p = cispi(phi)
    map(x -> p * x, Matrix2x2(su2))
end

struct SU2ParamScaled{T}
    uabs2::T
    alpha_u_scaled::T
    alpha_t_scaled::T
end

struct SU2Param1{T}
    u::Complex{T}
    t::Complex{T}
end

struct SU2Param2{T}
    uabs::T
    alpha_u::T
    alpha_t::T
end

struct SU2Param3{T}
    alpha_u::T
    alpha_t::T
    gamma::T
end

function unitary_decompose(U::Matrix2x2, ::Type{T}) where {T<:SU2Param1}
    T(U[1], U[2])
end

function unitary_decompose(U::Matrix2x2, ::Type{T}) where {T<:SU2Param3}
    (; u, t) = unitary_decompose(U, SU2Param1)
    (gamma, alpha_u, alpha_t) = _unitary_decompose_special(u, t)
    T(alpha_u, alpha_t, gamma)
end

end # module Matrices2x2
