module Gates

using ..Matrices2x2: Matrix2x2
using ..Common: canonical
using ..CyclotomicRings: Domega
using ..Singletons: One, Imag, RootImag, InvRootTwo
using ..Utils: PRETTY

export benchmark_compose

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

# Doing this has no effect on performance or resource use.
# It is only for cleaner looking code.
# Hmm. now that I use Gate1(:H), is this necessary? or useful?
# const Hgate = Gate1{:H}()
# const Sgate = Gate1{:S}()
# const Tgate = Gate1{:T}()
# const Wgate = Gate1{:W}()
# const Xgate = Gate1{:X}()
# const Ygate = Gate1{:Y}()
# const Zgate = Gate1{:Z}()

"""
    Matrix2x2{T}(g::Gate1{Name}) where {T, Name}
    Matrix2x2{T}(::Type{Gate1{Name}}) where {T, Name}

These methods instantiate a `Matrix2x2{T}` representing thate gate `Name`.

If a method for multiplying `Gate1(Name)` by a `Matrix2x2` is not defined, a `MethodError` will be thrown.
"""
Matrix2x2{T}(::Type{Gate1{Name}}) where {T, Name} = Matrix2x2{T}((Gate1{Name}() * one(Matrix2x2{T}))...,)
Matrix2x2{T}(g::Gate1{Name}) where {T, Name} = Matrix2x2{T}(typeof(g))

"""
    compose(gates::AbstractString; chunklen=300)

Compute composition of the gates in `gates`.

Runs of gates of length `chunklen` are computed using `Int64` as the primitive type.
The results are then composed with primitive type `BigInt`.

This function is not type-stable. If the length of `gates` is less than or equal to
`chunklen`, then no conversion to `BigInt` takes place.
# Examples
```jldoctest
julia> compose("TSHTHTHTHT")
2×2 Matrix2x2{Domega{Int}}:
1/2² ω³ + -1/2² ω² + 3/2² ω + 1/2² ω⁰   -1/2² ω³ + 1/2² ω² + 1/2² ω + 1/2² ω⁰
1/2² ω³ + 1/2² ω² + 1/2² ω + -1/2² ω⁰  -1/2² ω³ + -3/2² ω² + 1/2² ω + -1/2² ω⁰
```
"""
function compose(gates::AbstractString; chunklen=300)
    gates = reverse(gates)
    length(codeunits(gates)) <= chunklen && return compose_one(gates, false)
    chunks = reverse(chunkstring(gates, chunklen))
    mats = [map(Domega{BigInt}, compose_one(chunk, false)) for chunk in chunks]
    canonical(prod(mats))
end

"""
    chunkstring(s, chunklen::Integer)

Return an array of chunks of string `s`, each of length `chunklen`,
except the last one, which may be shorter. Satisfies
```
join(chunkstring(s)) == s
```
"""
function chunkstring(s, chunklen::Integer)
    strs = String[]
    i = 1
    while true
        cend = i + chunklen - 1
        cend = cend < length(s) ? cend : length(s)
        push!(strs, s[i:cend])
        cend == length(s) && break
        i = cend + 1
    end
    return strs
end

"""
    compose_one(gates::AbstractString, rev::Bool=true; reduce_fractions=true)
    compose_one(gates::AbstractString; reduce_fractions=true)

Compose gates in string `gates`. This is meant to compute the composition for a chunk
of a longer string of gates. The chunking and recombining is done by `compose`.
The reason we compute by chunks is that the chunks are small enough that computation
can be done with 64-bit (i.e. fast) integers. Then the resulting matrices are converted to
`BigInt`, and a final composition of these matrices is peformed.

Compute composition of the gates in `gates`. If `rev` is `true`, the composition will be from *right to
left*.

Following is disabled:
`map_func` is a map from `Symbol`s to matrices. If `reduce_fractions` is `true`
then reduce fractions in `Dyadic`s after each matrix multiplication.

`reduce_fractions` reduces the maximum values of intermediate numbers allowing computation
of longer compositions with smaller data types.
"""
function compose_one(gates::AbstractString, rev::Bool=true; reduce_fractions=true)
    gates = rev ? reverse(gates) : gates
    result = one(Matrix2x2{Domega{Int}})
    reduce_func = reduce_fractions ? canonical : identity
    for gate in codeunits(gates)
        new_result = _apply_gate(Symbol(Char(gate)), result)
        isnothing(new_result) && error(lazy"unknown gate $(Symbol(Char(gate)))")
        result = reduce_func(new_result)
    end
    return result
end

# Not used at the moment
# # `map_func` allows supplying gates that are not "built in"
# function compose_one(gates::AbstractString, map_func; reduce_fractions=true)
#     result = one(Matrix2x2{Domega{Int}})
#     reduce_func = reduce_fractions ? canonical : identity
#     for gate in codeunits(gates)
#         sgate = Symbol(Char(gate))
#         new_result = _apply_gate(sgate, result)
#         if isnothing(new_result)
#             new_result = map_func(sgate) * result
#         end
#         result = reduce_func(new_result)
#     end
#     return result
# end

# TODO: There are some packages that do this spliting for you.
# But I don't think general enough.
@inline function _apply_gate(gate::Symbol, matrix)
    if gate === :H
        Gate1(:H) * matrix
    elseif gate === :S
        Gate1(:S) * matrix
    elseif gate === :T
        Gate1(:T) * matrix
    elseif gate === :X
        Gate1(:X) * matrix
    elseif gate === :W
        Gate1(:W) * matrix
    else
        nothing
    end
end

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

"""
    count_gates(gates::AbstractString)::Dict{Char,Int}

Return a count `Dict` counting the number of each distinct gate in `gates`.

More precisely, this returns counts for each distinct character in `gates`.

# Example
```jldoctest
julia> count_gates("TSHTHTHTHT")
Dict{Char, Int64} with 3 entries:
  'H' => 4
  'T' => 5
  'S' => 1
```
"""
function count_gates(gates::AbstractString)
    counts = Dict{Char, Int}()
    for g in gates
        c = get(counts, g, 0)
        counts[g] = c + 1
    end
    return counts
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

# This mistake is very easy to make, causes bugs.
function Base.:*(::Type{T}, m::Matrix2x2) where {T <: Gate1{V}} where V
    throw(ArgumentError(lazy"Attempted matrix multiplication with a matrix type $(T), not a matrix value $(T)()"))
end

"""
    get_theta(m::Matrix2x2)

Find `theta` from a Z-rotation matrix `m` with possible global phase.

Assume `m` is diagonal, with `m[1] = cis(-theta/2 + phi)`
and `m[4] = cis(theta/2 + phi)`. Return `theta`.
"""
get_theta(m::Matrix2x2) = angle(m[4] / m[1])

get_theta(m::Matrix2x2{<:Domega}) = get_theta(big(m))

"""
    get_global_phase(m::Matrix2x2)

Find global phase of a Z-rotation matrix `m` with possible global phase.

Assume `m` is diagonal, with `m[1] = cis(-theta/2 + alpha)` and `m[4] = cis(theta/2 + alpha)`. Return `alpha`.

It's not possible to distinguish global phases that differ by addition of π. This is because
2α is extracted, and 2(α + π) will give the same argument of the phase factor.
"""
get_global_phase(m::Matrix2x2) = angle(m[4] * m[1]) / 2
get_global_phase(m::Matrix2x2{<:Domega}) = get_global_phase(big(m))

"""
    correct_global_phase(m::Matrix2x2)

This may be correct, or off by a global factor of `-1`.
"""
function correct_global_phase(m::Matrix2x2)
    alpha = get_global_phase(m)
    cis(-alpha) * big(m)
end

using LinearAlgebra: opnorm

"""
    rotation_error(m::matrix2x2, alpha::number; normf=opnorm)

compute the error in the approximate z-rotation gate `m` using the norm function `normf`.

`alpha` is the expected rotation angle. The RZ gate is `diagm([cis(-alpha/2), cis(alpha/2)])`
A global phase on the input matrix is corrected.
"""
function rotation_error(m::Matrix2x2, alpha::Number; normf=opnorm)
    expected_m = RZ(alpha)
    mb = big(m)
    mbc = correct_global_phase(mb)
    mdiff = mbc - expected_m
    msum = mbc + expected_m
    # We can only determine the global phase up to a term π.
    # So we try both.
    return min(normf(mdiff), normf(msum))
end

# using ..Dyadics: Dyadic
# function TSH()
#     T = Dyadic{Int64, Int64}
#     DT = Domega{Int64}
#     ab = DT((T(0, 0), T(1, 1), T(0, 0), T(-1, 1)))
#     c = DT((T(-1, 1), T(0, 0), T(1, 1), T(0, 0)))
#     d = DT((T(1, 1), T(0, 0), T(-1, 1), T(0, 0)))
#     Matrix2x2(ab, ab, c, d)
# end

end # module Gates
