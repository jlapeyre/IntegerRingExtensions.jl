module Compose

using DispatchDoctor: @stable, @unstable

using ..Matrices2x2: Matrix2x2, GPID, ScaleMatrix2x2
import ..Matrices2x2: get_theta
import ..RingMatrices: scalematrix
import ..Gates: Gate1, RZ
using ..Common: canonical
using ..CyclotomicRings: DOmega, ZOmega
using ..Utils: PRETTY, pretty
using ..Angles: random_angle

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
2×2 Matrix2x2{DOmega{Int}}:
1/2² ω³ + -1/2² ω² + 3/2² ω + 1/2² ω⁰   -1/2² ω³ + 1/2² ω² + 1/2² ω + 1/2² ω⁰
1/2² ω³ + 1/2² ω² + 1/2² ω + -1/2² ω⁰  -1/2² ω³ + -3/2² ω² + 1/2² ω + -1/2² ω⁰
```
"""
function compose(gates::AbstractString; chunklen=300, reduce_fractions=true)
    gates = reverse(gates)
    if length(codeunits(gates)) <= chunklen
        return scalematrix(compose_one(gates, false; reduce_fractions=reduce_fractions))
    end
    chunks = reverse(chunkstring(gates, chunklen))
    mats = [map(DOmega{BigInt}, compose_one(chunk, false)) for chunk in chunks]
    res = canonical(prod(mats))
    scalematrix(res)
end

# Matrix is numerically equal to `compose`, but differs sometimes by moving one
# factor of sqrt(2) from coeff to matrix.
# compose_scale is also much slower. I implemented in the hope that it would
# be faster
function compose_scale(gates::AbstractString; chunklen=300, reduce_fractions=true)
    gates = reverse(gates)
    if length(codeunits(gates)) <= chunklen
        return compose_one_scale(gates, false; reduce_fractions=reduce_fractions)
    end
    chunks = reverse(chunkstring(gates, chunklen))
    mats = [map(ZOmega{BigInt}, compose_one_scale(chunk, false)) for chunk in chunks]
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
# Using  scalematrix is a huge PITA
#    matrix = scalematrix(one(Matrix2x2{DOmega{Int}}))
    matrix = one(Matrix2x2{DOmega{Int}})
    reduce_func = reduce_fractions ? canonical : identity
    for gate in codeunits(gates)
        new_matrix = _apply_gate(Symbol(Char(gate)), matrix)
        isnothing(new_matrix) && error(lazy"unknown gate $(Symbol(Char(gate)))")
        matrix = reduce_func(new_matrix)
    end
     return matrix
end

# Using  scalematrix is a huge PITA
# For now, it's a separate function
function compose_one_scale(gates::AbstractString, rev::Bool=true; reduce_fractions=true)
    gates = rev ? reverse(gates) : gates
    matrix = scalematrix(one(Matrix2x2{DOmega{Int}}))
    reduce_func = reduce_fractions ? canonical : identity
    for gate in codeunits(gates)
        new_matrix = _apply_gate(Symbol(Char(gate)), matrix)
        isnothing(new_matrix) && error(lazy"unknown gate $(Symbol(Char(gate)))")
        matrix = reduce_func(new_matrix)
    end
    return matrix
end

function _alt_apply_gate(gates::AbstractString, rev::Bool=true; reduce_fractions=true)
    gates = rev ? reverse(gates) : gates
    matrix = one(Matrix2x2{DOmega{Int}})
    reduce_func = reduce_fractions ? canonical : identity
    for gate in codeunits(gates)
        gate = Symbol(Char(gate))
        new_result =
            if gate === :H
                Gate1{:H}() * matrix
            elseif gate === :T
                Gate1{:T}() * matrix
            elseif gate === :S
                Gate1{:S}() * matrix
            elseif gate === :Q  # We restrict gates in string to one char. So Q == Tdg
                Gate1{:Tdg}() * matrix
            elseif gate === :X
                Gate1{:X}() * matrix
            elseif gate === :W
                Gate1{:W}() * matrix
            elseif gate === :I
                Gate1{:I}() * matrix
            elseif gate === :Y
                Gate1{:Y}() * matrix
            elseif gate === :Z
                Gate1{:Z}() * matrix
            else
                throw(ArgumentError(lazy"Unsupported gate $gate"))
#                nothing
            end
 #       isnothing(new_result) && error(lazy"unknown gate $gate")
        matrix = reduce_func(new_result)
    end
    return matrix
end

# Only for testing.
# This incurs no allocation. OTOH, compose("HT"^n) incurs allocation for each iteration.
# I don't understand exactly why. But it must be due to the conditional
# An optimization would be to split the input string into sequences, like
# "HT"^n, etc. and apply these with a function like this one.
function apply_HT(n)
    matrix = one(Matrix2x2{DOmega{Int}})
    reduce = canonical
    for _ in 1:n
        matrix = reduce(Gate1{:H}() * matrix)
        matrix = reduce(Gate1{:T}() * matrix)
    end
    matrix
end

# Split using conditionals to avoid dynamic dispatch.
# i.e. don't dispatch on type of gate at runtime.
@inline function _apply_gate(gate::Symbol, matrix)
    if gate === :H
        Gate1{:H}() * matrix
    elseif gate === :T
        Gate1{:T}() * matrix
    elseif gate === :S
        Gate1{:S}() * matrix
    elseif gate === :Q  # We restrict gates in string to one char. So Q == Tdg
        Gate1{:Tdg}() * matrix
    elseif gate === :X
        Gate1{:X}() * matrix
    elseif gate === :W
        Gate1{:W}() * matrix
    elseif gate === :I
        Gate1{:I}() * matrix
    elseif gate === :Y
        Gate1{:Y}() * matrix
    elseif gate === :Z
        Gate1{:Z}() * matrix
    else
        throw(ArgumentError(lazy"Unsupported gate $gate"))
    end
end

#     if gate === :H
#         Gate1(:H) * matrix
#     elseif gate === :T
#         Gate1(:T) * matrix
#     elseif gate === :S
#         Gate1(:S) * matrix
#     elseif gate === :Q  # We restrict gates in string to one char. So Q == Tdg
#         Gate1(:Tdg) * matrix
#     elseif gate === :X
#         Gate1(:X) * matrix
#     elseif gate === :W
#         Gate1(:W) * matrix
#     elseif gate === :I
#         Gate1(:I) * matrix
#     elseif gate === :Y
#         Gate1(:Y) * matrix
#     elseif gate === :Z
#         Gate1(:Z) * matrix
#     else
#         throw(ArgumentError(lazy"Unsupported gate $gate"))
# #        matrix
# #        nothing
#     end
# end

"""
    get_theta(m::Matrix2x2)

Find `theta` from a Z-rotation matrix `m` with possible global phase.

Assume `m` is diagonal, with `m[1] = cis(-theta/2 + phi)`
and `m[4] = cis(theta/2 + phi)`. Return `theta`.
"""
function get_theta(m::Matrix2x2{Complex{T}}) where {T <:AbstractFloat}
    zsq = m[4] / m[1]
    angle(zsq)
end

get_theta(m::Matrix2x2) = get_theta(big(m))
get_theta(m::ScaleMatrix2x2) = get_theta(big(m))

"""
    get_global_phase(m::Matrix2x2)

Find global phase of a Z-rotation matrix `m` with possible global phase.

Assume `m` is diagonal, with `m[1] = cis(-theta/2 + alpha)` and `m[4] = cis(theta/2 + alpha)`. Return `alpha`.

It's not possible to distinguish global phases that differ by addition of π. This is because
2α is extracted, and 2(α + π) will give the same argument of the phase factor.
"""
get_global_phase(m::Matrix2x2{Complex{T}}) where {T <:AbstractFloat} = angle(m[4] * m[1]) / 2
get_global_phase(m::Matrix2x2) = get_global_phase(big(m))
get_global_phase(m::ScaleMatrix2x2) = get_global_phase(big(m))

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

function rotation_error_GPID(m::Matrix2x2, alpha::Number)
    expected_m = RZ(alpha)
    mb = big(m)
    return GPID(mb, expected_m)
end

# # I do not have a good name for this!
# """
#     Uapprox(theta, alpha, beta)::Matrix2x2

# Return a one-qubit unitary of the following type
# ```
# u    -t^*

# t    u^*
# ```
# where `u = cis(alpha)cos(theta)` and `t = cis(beta)sin(theta)`

# A more restricted case of this class of unitaries are those obtainable with
# a finite string of `H`, `S` and `cispi(1/4)`.
# """
# function Uapprox(theta, alpha, beta)
#     u = cis(alpha) * cos(theta)
#     t = cis(beta) * sin(theta)
#     Matrix2x2(u, t, -conj(t), conj(u))
# end

# function random_Uapprox(::Type{T}=Float64) where {T}
#     (theta, alpha, beta) = random_angle(T, 3)
#     Uapprox(theta, alpha, beta)
# end

# function isUapprox(m::Matrix2x2; kwargs...)
#     (u, t, mtconj, uconj) = m.data
#     isapprox(conj(u), uconj; kwargs...) || return false
#     isapprox(-conj(t), mtconj; kwargs...) || return false
#     isapprox(abs2(u) + abs2(t), 1; kwargs...) || return false
#     return true
# end

# using ..Dyadics: Dyadic
# function TSH()
#     T = Dyadic{Int64, Int64}
#     DT = DOmega{Int64}
#     ab = DT((T(0, 0), T(1, 1), T(0, 0), T(-1, 1)))
#     c = DT((T(-1, 1), T(0, 0), T(1, 1), T(0, 0)))
#     d = DT((T(1, 1), T(0, 0), T(-1, 1), T(0, 0)))
#     Matrix2x2(ab, ab, c, d)
# end

end # module Compose
