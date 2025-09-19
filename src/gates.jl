module Gates

using DispatchDoctor: @unstable

using ..Utils: PRETTY
import ..Matrices2x2: Matrix2x2, SU2, ScaleMatrix2x2
using ..RootOnes: omega, Omega
using ..CyclotomicRings: coeffs, div_half, mul_root_two, ZOmega

#using ..Singletons: Imag, RootImag, InvRootTwo, InvTwo, RootTwo
using ..Singletons: InvRootTwo, InvTwo, RootTwo

const RootImag = Omega(1)
const Imag = Omega(2)

export Gate1

const _onez = one(ZOmega)
const I2x2 = Matrix2x2(_onez, 0, 0, _onez)

"""
    X::Matrix2x2{ZOmega{Int64}}

The Pauli `X` matrix.
```
2×2 Matrix2x2{ZOmega{Int64}}:
 0  ω⁰
ω⁰   0
```
"""
const X = Matrix2x2(0, _onez, _onez, 0)

const Y = Matrix2x2(0, Imag * _onez, -(Imag * _onez), 0)
const Z = Matrix2x2(_onez, 0, 0, -_onez)
const S = Matrix2x2(_onez, 0, 0, omega^2)
const T = Matrix2x2(_onez, 0, 0, omega)
const H = InvRootTwo * (X + Z)

"""
    SX::Matrix2x2{DOmega{Int64}}

The square root of the Pauli `X` matrix.
```
2×2 Matrix2x2{DOmega{Int64}}:
 -1/2ω⁰ + 1/2ω²  -1/2ω⁰ + -1/2ω²
-1/2ω⁰ + -1/2ω²   -1/2ω⁰ + 1/2ω²
```
"""
const SX = Imag * (InvRootTwo * (omega * I2x2) + InvRootTwo * (omega^3 * X))
const SY = InvRootTwo * (omega * Matrix2x2(_onez, _onez, -_onez, _onez))

"""
    SH::Matrix2x2{DOmega{Int64}}

The square root of the `2x2` Hadamard matrix
```
2×2 Matrix2x2{DOmega{Int64}}:
1/2ω⁰ + 1/2ω² + -1/2ω³          -1/2ω³
        -1/2ω³          1/2ω⁰ + 1/2ω² + 1/2ω³
```
"""
const SH = InvTwo * Matrix2x2(-ZOmega(omega^3) + (RootTwo * ZOmega(omega)), -ZOmega(omega^3), -ZOmega(omega^3), ZOmega(omega^3) + (RootTwo * ZOmega(omega)))

commutator(x, y) = x * y - y * x

"""
    struct Gate1{Name} end

Represents a one-qubit gate with name `Name`

Be careful using these. It is not hard to end up with run-time dispatch, and a concomitant
huge performance hit.

# Examples
```jldoctest
julia> map(name -> Gate1(name), (:H, :S, :T, :X, :W))
(Gate1(:H), Gate1(:S), Gate1(:T), Gate1(:X), Gate1(:W))

julia> Gate1(:H) * one(Matrix2x2{DOmega{Int}})
2×2 Matrix2x2{DOmega{Int64}}:
1/2 ω + -1/2 ω³  1/2 ω + -1/2 ω³
1/2 ω + -1/2 ω³  -1/2 ω + 1/2 ω³

julia> Matrix2x2{DOmega{Int}}(Gate1(:H))
2×2 Matrix2x2{DOmega{Int64}}:
1/2 ω + -1/2 ω³  1/2 ω + -1/2 ω³
1/2 ω + -1/2 ω³  -1/2 ω + 1/2 ω³

julia> Matrix2x2{Float64}(Gate1(:H))
2×2 Matrix2x2{Float64}:
0.707107   0.707107
0.707107  -0.707107

julia> Matrix2x2{DOmega{Int}}(Gate1(:T))
2×2 Matrix2x2{DOmega{Int64}}:
 ω⁰   0
 0    ω
```
"""
struct Gate1{Name}
    @unstable function Gate1(name::Symbol)
        new{name}()
    end
    function Gate1{name}() where {name}
        isa(name, Symbol) || throw(MethodError(Gate1, name))
        new{name}()
    end
    function Gate1(name::Val{S}) where S
        new{S}()
    end
end

#Gate1(name::Symbol) = Gate1{name}()

function Base.show(io::IO, ::PRETTY, ::Gate1{Name}) where {Name}
    print(io, "Gate1{$(repr(Name))}()")
end
Base.show(io::IO, g::Gate1) = show(io, PRETTY(), g)

"""
    Matrix2x2{T}(g::Gate1{Name}) where {T, Name}
    Matrix2x2{T}(::Type{Gate1{Name}}) where {T, Name}

These methods instantiate a `Matrix2x2{T}` representing thate gate `Name`.

If a method for multiplying `Gate1(Name)` by a `Matrix2x2` is not defined, a `MethodError` will be thrown.
"""
Matrix2x2{T}(::Type{Gate1{Name}}) where {T, Name} = Matrix2x2{T}((Gate1{Name}() * one(Matrix2x2{T}))...,)
Matrix2x2{T}(g::Gate1{Name}) where {T, Name} = Matrix2x2{T}(typeof(g))


"""
    RZ(theta)

This is the Z rotation matrix that is synthesized by `gridsynth`
 as well as algorithms in other papers.
"""
function RZ(theta)
    z = zero(theta)
    t2 = theta/2
    Matrix2x2(cis(-t2), z, z, cis(t2))
end

# function random_RZ(::Type{T}=Float64) where {T}
#     theta = random_angle(T)
#     return RZ(theta)
# end

Base.:*(::Gate1{:W}, m::Matrix2x2) =  map(x-> RootImag * x, m)
Base.:*(m::Matrix2x2, ::Gate1{:W}) = Gate1{:W}() * m

function Base.:*(::Gate1{:X}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(b, a, d, c)
end

function Base.:*(::Gate1{:Y}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(- Imag *b, Imag * a, -(Imag * d), Imag * c)
end

function Base.:*(::Gate1{:Z}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, -b, c, -d)
end

function Base.:*(m::Matrix2x2, ::Gate1{:X})
    (a,b,c,d) = m.data
    Matrix2x2(c, d, a, b)
end

function Base.:*(::Gate1{:S}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, Imag * b, c, Imag * d)
end

function Base.:*(::Gate1{:Sdg}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, -(Imag * b), c, -(Imag * d))
end

function Base.:*(::Gate1{:H}, m::Matrix2x2)
    (a,b,c,d) = m.data
    s = InvRootTwo
    Matrix2x2(s*(a + b), s*(a-b), s*(c+d), s*(c-d))
end

# We could also use `omega` instead of `RootImag`.
function Base.:*(::Gate1{:T}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, RootImag * b, c, RootImag * d)
end

function Base.:*(::Gate1{:Tdg}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, inv(omega) * b, c, inv(omega) * d)
end

function Base.:*(::Gate1{:I}, m::Matrix2x2)
    m
end

#################

Base.:*(::Gate1{:W}, m::ScaleMatrix2x2) =  ScaleMatrix2x2(omega * m.m, m.s)
Base.:*(m::ScaleMatrix2x2, ::Gate1{:W}) = Gate1{:W}() * m

# Trying different, equivalent ways to do this multiplication.
# This seems to be correct. Agrees with version using Matrix2x2.
# However the coefficients are 10^14 times larger. This is unacceptable and uneeded
function Base.:*(::Gate1{:H}, m::ScaleMatrix2x2)
    (a,b,c,d) = m.m.data

    x = a + b
    y = c + d
    if all(iseven, coeffs(x)) && all(iseven, coeffs(y))
#        print("1")
        (a2, b2, c2, d2) = (a + b, a - b, c + d, c - d)
        (an, bn, cn, dn) = map(x -> mul_root_two(div_half(x)), (a2, b2, c2, d2))
        return ScaleMatrix2x2(Matrix2x2(an, bn, cn, dn), m.s)
    end

    s_new = InvRootTwo * m.s
    m_new = Matrix2x2(x, a - b, y, c - d)
    ScaleMatrix2x2(m_new, s_new)
end

function Base.:*(g::Gate1, m::ScaleMatrix2x2)
    ScaleMatrix2x2(typeof(g)() * m.m, m.s)
end

# function Base.:*(::Gate1{:X}, m::ScaleMatrix2x2)
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, (b, a, d, c))
# end

# function Base.:*(::Gate1{:Y}, m::ScaleMatrix2x2)
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, Matrix2x2(b, Imag * a, -(Imag * d), c))
# end

# function Base.:*(::Gate1{:Z}, m::ScaleMatrix2x2)
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, Matrix2x2(a, -b, c, -d))
# end

# function Base.:*(m::ScaleMatrix2x2, ::Gate1{:X})
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, Matrix2x2(c, d, a, b))
# end

# function Base.:*(::Gate1{:S}, m::ScaleMatrix2x2)
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, Matrix2x2(a, Imag * b, c, Imag * d))
# end

# function Base.:*(::Gate1{:Sdg}, m::ScaleMatrix2x2)
#     (a,b,c,d) = m.m.data
#     ScaleMatrix2x2(m.s, Matrix2x2(a, -(Imag * b), c, -(Imag * d)))
# end

# function Base.:*(::Gate1{:H}, m::ScaleMatrix2x2)
#     news = InvRootTwo * m.s
#     ScaleMatrix2x2(Gate1(:H) * m.m, news)
# end

# # We could also use `omega` instead of `RootImag`.
# function Base.:*(::Gate1{:T}, m::ScaleMatrix2x2)
#     ScaleMatrix2x2(m.s, Gate1(:T) * m.m)
# end

# function Base.:*(g::Gate1, m::ScaleMatrix2x2)
#     ScaleMatrix2x2(typeof(g)() * m.m, m.s)
# end

# function Base.:*(::Gate1{:Tdg}, m::Matrix2x2)
#     ScaleMatrix2x2(m.s, Gate1(:Tdg) * m.m)
# end

# function Base.:*(::Gate1{:I}, m::Matrix2x2)
#     ScaleMatrix2x2(m.s, Gate1(:I) * m.m)
# end

#################


# This mistake is very easy to make, causes bugs.
function Base.:*(::Type{T}, m::Matrix2x2) where {T <: Gate1{V}} where V
    throw(ArgumentError(lazy"Attempted matrix multiplication with a matrix type $(T), not a matrix value $(T)()"))
end

function SU2(::Type{Gate1{:V1}}, ::Type{T}=Float64) where {T}
    d = T(1)/sqrt(T(5))
    SU2(complex(d), Complex(0, T(2)) * d)
end

function SU2(::Type{Gate1{:V2}}, ::Type{T}=Float64) where {T}
    d = T(1)/sqrt(T(5))
    SU2(complex(d), complex(0, T(2)) * d)
end

function SU2(::Type{Gate1{:V3}}, ::Type{T}=Float64) where {T}
    d = T(1)/sqrt(T(5))
    SU2(Complex(d, 2*d), complex(zero(d)))
end

function Matrix2x2(::Type{Gate1{:H}}, ::Type{T}=Float64) where {T}
    d = T(1)/sqrt(T(2))
    Matrix2x2(d, d, d, -d)
end

function Matrix2x2(::Type{Gate1{:S}}, ::Type{T}=Float64) where {T}
    Matrix2x2(complex(one(T)), zero(T), zero(T), Complex(zero(T), one(T)))
end


end # module Gates
