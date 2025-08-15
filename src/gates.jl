module Gates

using ..Utils: PRETTY
using ..Matrices2x2: Matrix2x2
using ..Singletons: Imag, RootImag, InvRootTwo

export Gate1

"""
    struct Gate1{Name} end

Represents a one-qubit gate with name `Name`

Be careful using these. It is not hard to end up with run-time dispatch, and a concomitant
huge performance hit.

# Examples
```jldoctest
julia> map(name -> Gate1(name), (:H, :S, :T, :X, :W))
(Gate1(:H), Gate1(:S), Gate1(:T), Gate1(:X), Gate1(:W))

julia> Gate1(:H) * one(Matrix2x2{Domega{Int}})
2×2 Matrix2x2{Domega{Int64}}:
1/2 ω + -1/2 ω³  1/2 ω + -1/2 ω³
1/2 ω + -1/2 ω³  -1/2 ω + 1/2 ω³

julia> Matrix2x2{Domega{Int}}(Gate1(:H))
2×2 Matrix2x2{Domega{Int64}}:
1/2 ω + -1/2 ω³  1/2 ω + -1/2 ω³
1/2 ω + -1/2 ω³  -1/2 ω + 1/2 ω³

julia> Matrix2x2{Float64}(Gate1(:H))
2×2 Matrix2x2{Float64}:
0.707107   0.707107
0.707107  -0.707107

julia> Matrix2x2{Domega{Int}}(Gate1(:T))
2×2 Matrix2x2{Domega{Int64}}:
 ω⁰   0
 0    ω
```
"""
struct Gate1{Name} end

Gate1(name) = Gate1{name}()
function Base.show(io::IO, ::PRETTY, ::Gate1{Name}) where {Name}
    print(io, "Gate1($(repr(Name)))")
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

function random_RZ(::Type{T}=Float64) where {T}
    theta = random_angle(T)
    return RZ(theta)
end

Base.:*(::Gate1{:W}, m::Matrix2x2) =  map(x-> RootImag * x, m)
Base.:*(m::Matrix2x2, ::Gate1{:W}) = Gate1(:W) * m

function Base.:*(::Gate1{:X}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(b, a, d, c)
end

function Base.:*(::Gate1{:Y}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(b, Imag * a, -(Imag * d), c)
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

function Base.:*(::Gate1{:H}, m::Matrix2x2)
    (a,b,c,d) = m.data
    s = InvRootTwo
    Matrix2x2(s*(a + b), s*(a-b), s*(c+d), s*(c-d))
end

function Base.:*(::Gate1{:T}, m::Matrix2x2)
    (a,b,c,d) = m.data
    Matrix2x2(a, RootImag * b, c, RootImag * d)
end

function Base.:*(::Gate1{:I}, m::Matrix2x2)
    m
end

# This mistake is very easy to make, causes bugs.
function Base.:*(::Type{T}, m::Matrix2x2) where {T <: Gate1{V}} where V
    throw(ArgumentError(lazy"Attempted matrix multiplication with a matrix type $(T), not a matrix value $(T)()"))
end


end # module Gates
