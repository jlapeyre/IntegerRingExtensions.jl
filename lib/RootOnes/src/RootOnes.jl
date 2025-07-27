"""
module RootOnes

- [`struct RootOne`](@ref "RootOne") – represents roots of unity
- [`const Omega`](@ref "Omega") – Alias for `RootOne{8}`, the eighth roots of unity.
- [`const omega`](@ref "omega") – `Omega(1)`, the principal eighth roots of unity.

The following methods from `RingExtensionsCommon` are implemented for `RootOne`.

[`imaginary`](@ref), [`sqrt_imaginary`](@ref), [`isunit`](@ref), [`conj_root_two`](@ref), [`norm_root_two`](@ref)

A few other methods
[`rand`](@ref)

Examples of several more methods are given in [`RootOne`](@ref).
"""
module RootOnes

# IDK how to use Revise with this. Always errors on every edit
# @stable module RootOnes

import Base: convert, show

import RingExtensionsUtils: subscript, superscript

# Items defined in this module
export RootOne, Omega, omega, isprimitive

# Re-export imported functions
import RingExtensionsCommon: imaginary, sqrt_imaginary, isunit, conj_root_two, norm_root_two
export imaginary, sqrt_imaginary, isunit, conj_root_two, norm_root_two

import Random: Random, rand
export rand

########################
####
#### RootOne
####
########################

# This struct implements the cyclic group of order N,
# plus a few methods particular to the interpretation as roots of unity.
"""
    struct RootOne{N}
    RootOne{N}(k::Integer)

`N`th roots of unity in the complex numbers.

`RootOne{N}(k)` represents `exp(2πi k / N)`.

`k` is stored as `mod(k, N)`, which takes values in `(0,...,N-1)`.

We call `exp(2πi / N)`, represented by `RootOne{N}(1)`, the *principal* root.

The following aliases are defined:
```julia
const Omega = RootOne{8}
const omega = Omega(1) # Principal eighth root of unity
```

See also [`Omega`](@ref) and [`omega`](@ref).

!!! warning
    `N` should be a literal or a `const`. Otherwise performance is severely degraded.

## Examples

These examples include methods that have no docstrings

The principal `N`ᵗʰ root is constructed like this
```jldoctest
julia> RootOne{8}(1)
ω₈

julia> RootOne{3}() # Omit the argument for principal root
ω₃
```

`k` is shifted by `N` such that `k ∈ (0, ..., N-1)`.
```jldoctest
julia> RootOne{8}(-1)
ω₈⁷

julia> omega^-1 === omega^7 === omega^15
true
```

```jldoctest
julia> one(Omega)
ω₈⁰

julia> one(RootOne{5})
ω₅⁰

julia> one(RootOne{8}) === one(RootOne{8}(5)) === RootOne{8}(0)
true

julia> RootOne(2, 6) == RootOne(1, 3)
true

julia> RootOne(1, 6) == RootOne(1, 3)
false
```

### Arithmetic

Unary minus is not a group operation.
```jldoctest
julia> -RootOne{8}(1) # Unary minus
ω₈⁵

julia> -RootOne{8}(5)
ω₈

julia> -RootOne{9}(1)
ERROR: ArgumentError: InexactError unary minus of type RootOne{9}
```

Multiplication and powers operate on elements of the complex numbers.
They may also be group operations.
```jldoctest
julia> RootOne(2, 8) * RootOne(3, 8) # Group operation
ω₈⁵

julia> RootOne{8}(1)^8
ω₈⁰

julia> isone(RootOne{8}(1)^8)
true

julia> omega^7 / omega
ω₈⁶
```

Multiplication between different groups is not a group operation and may return a value from a third group.
```jldoctest
julia> RootOne(1, 3) * RootOne(1, 2)
ω₆⁵
```

Some roots are real.
```jldoctest
julia> Float64(Omega(0))
1.0

julia> Float64(Omega(4))
-1.0
```

### Conversion
```jldoctest
julia> complex(omega)
0.7071067811865476 + 0.7071067811865475im

julia> float(omega)
0.7071067811865476 + 0.7071067811865475im

julia> complex(RootOne(1, 4))
0 + 1im

julia> float(RootOne(1, 4))
0.0 + 1.0im
```

### Functions on complex numbers
```jldoctest
julia> angle(omega)
0.7853981633974483

julia> angle(BigFloat, omega)
0.7853981633974483096156608458198757210492923498437764552437361480769541015715495

julia> log(omega)
0.0 + 0.7853981633974483im

julia> inv(omega)
ω₈⁷

julia> conj(omega) === inv(omega)
true

julia> abs(omega) === abs2(omega) === 1
true

julia> real(omega)
0.7071067811865476

julia> real(BigFloat, omega)
0.707106781186547524400844362104849039284835937688474036588339868995366239231051

julia> imag(omega^2)
1.0

julia> isreal(omega)
false

julia> isreal(omega^4)
true

julia> isreal(omega^0)
true

julia> isinteger(omega^0)
true
```
"""
struct RootOne{N} <: Number
    function RootOne{N}(k::Integer) where {N} # Annotate `Integer` for method ambiguities.
        N isa Integer || throw(ArgumentError(lazy"In RootOne{N}, N must be an integer"))
        N > 0 || throw(ArgumentError(lazy"In RootOne{N}, N must be > 0"))
        n = Int(N) # Not sure what this normalization is worth.
        new{n}(mod(k, n))
    end
    k::Int # Could use smaller `k`.
end

RootOne{N}() where {N} = RootOne{N}(1)

"""
    RootOne(k::Integer, N::Integer)

Return `RootOne{Int(N)}(k)`.

This is an alternative constructor.

### Examples
```jldoctest
julia> RootOne(3, 10)
ω₁₀³

julia> RootOne{10}(3)
ω₁₀³
```
"""
RootOne(k::Integer, N::Integer) = RootOne{Int(N)}(k)

"""
    const Omega = RootOne{8}

The eighth roots of unity.

### Examples

```julia-repl
julia> Omega()
ω₈

julia> Omega(1)
ω₈

julia> Omega(0)
ω₈⁰

julia> Omega(5)
ω₈⁵

julia> Omega(5) * Omega(3)
ω₈⁰
```
"""
const Omega = RootOne{8}

"""
    const omega = Omega(1)

The principal eighth root of unity `exp(2πi/8)`

### Examples
```jldoctest
julia> omega
ω₈

julia> omega^8
ω₈⁰

julia> (omega, omega^8)
(ω₈, ω₈⁰)

julia> isone(omega^8)
true
```
"""
const omega = Omega(1)

"""
    imaginary(::Type{RootOne{D}}) where {D}

The value of type `RootOne{D}` that represents the imaginary unit, 𝕚.

An error is thrown if no such value exists.

### Examples
```jldoctest
julia> imaginary(RootOne{4})
ω₄

julia> imaginary(RootOne{8})
ω₈²

julia> imaginary(RootOne{3})
ERROR: ArgumentError: Type RootOne{3} cannot represent 𝕚
```
"""
function imaginary(::Type{RootOne{D}}) where {D}
    (n, resid) = divrem(D, 4)
    resid == 0 || throw(ArgumentError(lazy"Type RootOne{$D} cannot represent 𝕚"))
    RootOne{D}(n)
end

"""
    sqrt_imaginary(::Type{RootOne{D}})::RootOne{D} where {D}

The value of type `RootOne{D}` that represents √𝕚, the principal square root of the imaginary unit.

An error is thrown if no such value exists.

### Examples
```jldoctest
julia> sqrt_imaginary(Omega)
ω₈

julia> sqrt_imaginary(RootOne{16})
ω₁₆²

julia> sqrt_imaginary(RootOne{12})
ERROR: ArgumentError: Type RootOne{12} cannot represent √𝕚
```
"""
function sqrt_imaginary(::Type{RootOne{D}}) where {D}
    (n, resid) = divrem(D, 8)
    resid == 0 || throw(ArgumentError(lazy"Type RootOne{$D} cannot represent √𝕚"))
    RootOne{D}(n)
end

function _print_omega(io::IO, N::Int, k::Int)
    print(io, "ω" * subscript(N))
    isone(k) || print(io, superscript(k))
    nothing
end

Base.show(io::IO, r::RootOne{N}) where {N} = _print_omega(io, N, r.k)
Base.show(io::IO, ::MIME"text/plain", r::RootOne{N}) where {N} = _print_omega(io, N, r.k)

Base.one(::Type{RootOne{N}}) where {N} = RootOne{N}(0)
Base.one(x::RootOne) = one(typeof(x))
Base.isone(x::RootOne) = iszero(x.k) # fallback is probably ok, too

# TODO: Consider letting Base throw method error. Or we throw it explicitly
function Base.zero(::Type{<:RootOne})
    throw(ArgumentError("Type RootOne has no additive identity"))
end
Base.zero(::T) where {T<:RootOne} = zero(RootOne)

RootOne{N}(r::RootOne{N}) where {N} = r

"""
    RootOne{M}(rt::RootOne{N}) where {M, N}

Convert `rt` to an instance of `RootOne{M}`.

An error is thrown if exact conversion is not possible

### Examples
```jldoctest
julia> RootOne{4}(RootOne{8}(6))
ω₄³

julia> RootOne{4}(RootOne{8}(1))
ERROR: ArgumentError: InexactError converting ω₈ to RootOne{4}
```
"""
function RootOne{M}(rt::RootOne{N}) where {M, N}
    if M < N
        (n, r) = divrem(N, M)
        r == 0 || throw(ArgumentError(lazy"InexactError converting $rt to RootOne{$M}"))
        (q, r) = divrem(rt.k, n)
        r == 0 || throw(ArgumentError(lazy"InexactError converting $rt to RootOne{$M}"))
        return RootOne{M}(q)
    else
        (n, r) = divrem(M, N)
        r == 0 || throw(ArgumentError(lazy"InexactError converting $rt to RootOne{$M}"))
        return RootOne{M}(rt.k * n)
    end
end

function Base.:(==)(r1::RootOne{N1}, r2::RootOne{N2}) where {N1,N2}
    r1.k // N1 == r2.k // N2
end

Base.hash(r::RootOne{N}, h::UInt) where {N} = hash(r.k // N, h)

function (::Type{T})(r::RootOne{M}; maybe::Bool=false) where {M, T <: Real}
    iszero(r.k) && return one(T)
    if !(iseven(M) && div(M, r.k) == 2)
        maybe && return nothing
        throw(ArgumentError(lazy"InexactError converting RootOne{$M}($r) to $T"))
    end
    return -one(T)
end

# TODO: How can I see this docstring?
"""
    Complex{T}(r::RootOne{M}; maybe=false) where {M, T <: Integer}

Return `r` as a `Complex{T<:Integer}`.

If no such conversion exists either throw an `ArgumentError`, if `maybe` is `false`,
or return `nothing`, if `maybe` is `true`.

### Examples
Some roots are Gaussian integers
```jldoctest
julia> Tuple((i, Complex{Int}(Omega(i))) for i in (0,2,4,6))
((0, 1 + 0im), (2, 0 + 1im), (4, -1 + 0im), (6, 0 - 1im))
```
"""
function (::Type{Complex{T}})(r::RootOne{M}; maybe::Bool=false) where {M, T <: Integer}
    CT = Complex{T}
    iszero(r.k) && return one(CT)
    (n, chk) = divrem(r.k << 2, M)
    if !iszero(chk)
        maybe && return nothing
        throw(ArgumentError(lazy"InexactError converting RootOne{$M}($r) to $T"))
    end
    (z, o) = (zero(T), one(T))
    n == 1 && return CT(z, o)
    n == 2 && return CT(-o, z)
    return CT(z, -o)
end

Base.float(r::RootOne) = float(complex(r))
Base.big(r::RootOne) = Complex{BigFloat}(r)

# We don't want these. Rather convert to: Complex{QuadraticRing{2, Dyadic{Int64, Int64}}}
# by defining method for `Base.complex` in quadratic_ring.jl.
# Base.complex(r::RootOne) = Complex(r)
# convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N} = Complex{T}(r)

# Base.Complex(r::RootOne) = Complex{Float64}(r)

Base.Complex{T}(r::RootOne{N}) where {T, N} = cispi((T(2) * T(r.k)) / T(N))

# Method for `Base.complex(r::RootOne{8})` defined in quadratic_ring.jl.
#Base.Complex(r::RootOne) = complex(r)
Base.Complex(r::RootOne) = Complex{Float64}(r)
Base.complex(r::RootOne) = Complex(r)

Base.angle(r::RootOne) = angle(Float64, r)
function Base.angle(::Type{T}, r::RootOne{N}) where {N, T}
    k = r.k <= N/2 ? r.k : r.k - N
    2 * T(pi) * T(k) / T(N)
end

Base.log(r::RootOne) = log(Float64, r)
function Base.log(::Type{T}, r::RootOne) where {T}
    Complex(zero(T), angle(T, r))
end

Base.adjoint(r::RootOne) = conj(r)
Base.conj(r::RootOne) = inv(r)
Base.abs2(r::RootOne) = 1
Base.abs(r::RootOne) = 1

# TODO: check conj_root_two(omega^2). Is this correct?
"""
    conj_root_two(r::RootOne{8})

Return `r` with √2 replaced by -√2.
"""
conj_root_two(r::RootOne{8}) = -r

"""
    norm_root_two(r::RootOne{8})

Return `r` time `conj_root_two(r)`.
"""
norm_root_two(r::RootOne{8}) = -r^2

Base.real(r::RootOne{N}) where {N} = real(Float64, r)
Base.real(::Type{T}, r::RootOne{N}) where {T, N} = cospi(T(2 * r.k) / N)
Base.imag(r::RootOne{N}) where {N} = imag(Float64, r)
Base.imag(::Type{T}, r::RootOne{N}) where {T, N} = sinpi(T(2 * r.k) / N)

function Base.isreal(r::RootOne{N}) where {N}
    r.k == 0 && return true  # r == 1
    iszero(mod(r.k << 1, N))  # r == -1
end

Base.isinteger(r::RootOne) = isreal(r)

#Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(mod(r1.k + r2.k, N))
Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k + r2.k)
Base.:/(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k - r2.k)

function Base.promote_rule(::Type{RootOne{N}}, ::Type{RootOne{M}}) where {N,M}
    return RootOne{lcm(N, M)}
end

# This is stable
function Base.:*(a::RootOne{N}, b::RootOne{M}) where {N,M}
    x, y = promote(a, b)
    return x * y
end

function Base.:/(r1::RootOne{N}, r2::RootOne{M}) where {N, M}
    r1 * inv(r2)
end

function Base.sqrt(r::RootOne{N}) where N
    k = r.k
    iseven(k) || throw(ArgumentError(lazy"Inexact error: sqrt($r)"))
    RootOne{N}(k >> 1)
end

#      Base.literal_pow(::typeof(^), r::RootOne{N}, ::Val{n}) where {N,n} = RootOne{N}(mod(r.k * mod(n, N), N))
Base.literal_pow(::typeof(Base.:^), r::RootOne{N}, ::Val{n}) where {n, N} = RootOne{N}(r.k * mod(n, N))
Base.:^(r::RootOne{N}, n::Integer) where {N} = RootOne{N}(r.k * mod(n, N))

# TODO: figure out how to attach and receive doc info correctly. (*sigh*)
# `? inv(::RootOne)`, etc. do not work as expected.
"""
    inv(r::RootOne{N})::RootOne{N}

The inverse of `r`.

`isone(r * inv(r))` is `true`.
"""
Base.inv(::RootOne) # Does not help
Base.inv(r::RootOne{N}) where {N} = RootOne{N}(N - r.k)

"""
    isunit(::RootOne)::Bool

Always returns `true`.

All elements of the group of `n`ᵗʰ roots of unity are invertible.
"""
isunit(::RootOne) = true

"""
    isprimitive(r::RootOne{N})::Bool where N

Return `true` if `r` is a primitive root.

`r` is a primitive root if and only if the first `N-1` powers of `r` are distinct.

### Examples

```jldoctest
julia> Tuple((k, isprimitive(omega^k)) for k in 1:8)
((1, true), (2, false), (3, true), (4, false), (5, true), (6, false), (7, true), (8, false))
```
"""
isprimitive(r::RootOne{N}) where N = isone(gcd(r.k, N))

# Have to be careful with call syntax
# Use this:
# Base.:-(RootOne{7}(3); maybe=true)
# Not this:
# -(RootOne{7}(3); maybe=true)
function Base.:-(r::RootOne{N}; maybe=false) where {N}
    if isodd(N)
        maybe && return nothing
        throw(ArgumentError(lazy"InexactError unary minus of type RootOne{$N}"))
    end
    RootOne{N}(r.k + (N >> 1))
end

# TODO: The represention as `Float64` is also exact. Which to choose
# TODO: Check why "duplicate" methods are present. Often for disambiguation.
# If this is so, comment these methods so they don't have to be checked later.
Base.real(::RootOne{1}) = 1
Base.imag(::RootOne{1}) = 0
Base.Complex(r::RootOne{1}) = Complex{Int}(r)
Base.Complex{T}(r::RootOne{1}) where {T<:Number} = Complex(one(T), zero(T))
Base.Complex{T}(r::RootOne{1}) where {T<:Integer} = Complex(one(T), zero(T))
Base.Complex{Int}(r::RootOne{1}) = Complex(1, 0)

Base.real(::Type{T}, r::RootOne{2}) where T <: Number = T(real(r))
Base.real(r::RootOne{2}) = Int(r)
Base.Complex(r::RootOne{2}) = Complex{Int}(r)

Base.Complex(r::RootOne{4}) = Complex{Int}(r)

"""
    Complex(r::RootOne{4})
    Complex{T}(r::RootOne{4}) where {T <: Integer}

Convert instances of `RootOne{4}` to `Complex{<:Integer}`

The first form returns `Complex{Int}`.

### Examples

```jldoctest
julia> [Complex(RootOne(i, 4)) for i in 0:3]
4-element Vector{Complex{Int64}}:
  1 + 0im
  0 + 1im
 -1 + 0im
  0 - 1im
```
"""
function Base.Complex{T}(r::RootOne{4}) where {T <: Integer}
    k = r.k
    o = one(T)
    z = zero(T)
    tup =
        if k == 0
            (o, z)
        elseif k == 1
            (z, o)
        elseif k == 2
            (-o, z)
        else k == 3
            (z, -o)
        end
    Complex(tup...)
end

"""
    rand([rng=default_rng()], RootOne{N}, [dims...]) where {N}

Return a random element or array of `RootOne{N}`.

### Examples
```julia-repl
julia> rand(Omega)
ω₈³

julia> rand(Omega, (2,2))
2×2 Matrix{Omega}:
 ω₈⁷  ω₈⁴
  ω₈  ω₈⁴
```
"""
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{RootOne{N}}) where {N}
    RootOne{N}(rand(rng, 0:N-1))
end

end # module RootOnes
