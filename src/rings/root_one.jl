module RootOnes

import ..Utils: subscript, superscript
import ..Common: sqrt_imaginary, imaginary
import Base: convert, show
import Random

export RootOne, Omega

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

`N`th roots of unity. `RootOne{N}(k)` is the `k`th power of the principal root.

`N` should always be a literal or a `const`. Otherwise performance is severely degraded.

`k` will be stored as `mod(k, N)`, which takes values in `(0,...,N-1)`.

The following alias is defined:
```julia
const Omega = RootOne{8}
```

# Examples

See also `Omega`.

```jldoctest

julia> RootOne{8}(1)
RootOne{8}(1)
```

```jldoctest
julia> RootOne{8}(-1)
RootOne{8}(7)

julia> -RootOne{8}(1) # Unary minus
RootOne{8}(5)

julia> -RootOne{8}(5)
RootOne{8}(1)

julia> -RootOne{9}(1)
ERROR: ArgumentError: InexactError unary minus of type RootOne{9}

julia> RootOne{8}(1)^8
RootOne{8}(0)

julia> isone(RootOne{8}(1)^8)
true

julia> one(RootOne{8})
RootOne{8}(0)

julia> one(RootOne{8}(5))
RootOne{8}(0)

julia> RootOne{8}(2) * RootOne{8}(3)
RootOne{8}(5)
```
"""
struct RootOne{N}
    function RootOne{N}(k) where {N}
        new{N}(mod(k, N))
    end
    k::Int # Could use smaller `k`.
end

"""
    RootOne{N}(n=1)

Return `RootOne{N}(1)`, the principal `n`th root of unity.
"""
RootOne{N}() where {N} = RootOne{N}(1)

"""
    Omega = RootOne{8}

The eighth roots of unity.

# Examples

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

julia> complex(Omega(3))
(0 + -1/2√2) + (0 + 1/2√2)im
```
"""
const Omega = RootOne{8}

"""
    imaginary(::Type{RootOne{D}}) where {D}

The value of type `RootOne{D}` that represents the imaginary unit.

An error is thrown if no such value exists.
"""
function imaginary(::Type{RootOne{D}}) where {D}
    (n, resid) = divrem(D, 4)
    resid == 0 || throw(ArgumentError(lazy"Type RootOne{$D} cannot represent √𝕚"))
    RootOne{D}(n)
end

"""
    sqrt_imaginary(::Type{RootOne{D}}) where {D}

The value of type `RootOne{D}` that represents the principal square root of the imaginary unit.

An error is thrown if no such value exists.
"""
function sqrt_imaginary(::Type{RootOne{D}}) where {D}
    (n, resid) = divrem(D, 8)
    resid == 0 || throw(ArgumentError(lazy"Type RootOne{$D} cannot represent √𝕚"))
    RootOne{D}(n)
end

function show(io::IO, ::MIME"text/plain", r::RootOne{N}) where {N}
    print(io, "ω" * subscript(N))
    isone(r.k) || print(io, superscript(r.k))
end

"""
    RootOne(k, N)

Construct `RootOne{N}(k)`.
"""
RootOne(k, N) = RootOne{N}(k)
Base.one(::Type{RootOne{N}}) where {N} = RootOne{N}(0)
Base.one(x::RootOne) = one(typeof(x))

function Base.:(==)(r1::RootOne{N1}, r2::RootOne{N2}) where {N1,N2}
    r1.k // N1 == r2.k // N2
end

(::Type{T})(r::RootOne{1}) where {T <: Number} = one(T)
(::Type{T})(r::RootOne{1}) where {T <: Complex} = one(T)
Base.Complex(r::RootOne{1}) = Complex{Float64}(r)

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
Base.complex(r::RootOne) = Complex(r)
Base.big(r::RootOne) = Complex{BigFloat}(r)
convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N} = Complex{T}(r)
Base.Complex{T}(r::RootOne{N}) where {T, N} = cispi((T(2) * T(r.k)) / T(N))
Base.Complex(r::RootOne) = Complex{Float64}(r)

Base.angle(r::RootOne) = angle(Float64, r)
function Base.angle(::Type{T}, r::RootOne{N}) where {N, T}
    k = r.k <= N/2 ? r.k : r.k - N
    2 * T(pi) * T(k) / T(N)
end

Base.adjoint(r::RootOne) = conj(r)
Base.conj(r::RootOne) = inv(r)
Base.abs2(r::RootOne) = 1
Base.abs(r::RootOne) = 1

Base.real(r::RootOne{N}) where {N} = real(Int, r)
Base.real(::Type{T}, r::RootOne{N}) where {T, N} = cospi(T(2 * r.k) / N)
Base.imag(r::RootOne{N}) where {N} = imag(Int, r)
Base.imag(::Type{T}, r::RootOne{N}) where {T, N} = sinpi(T(2 * r.k) / N)

function Base.isreal(r::RootOne{N}) where {N}
    r.k == 0 && return true  # r == 1
    iszero(mod(r.k << 1, N))  # r == -1
end

Base.isinteger(r::RootOne) = isreal(r)

Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k + r2.k)
Base.:/(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k - r2.k)

function Base.sqrt(r::RootOne{N}) where N
    k = r.k
    iseven(k) || throw(ArgumentError(lazy"Inexact error: sqrt($r)"))
    RootOne{N}(k >> 1)
end

Base.:^(r::RootOne{N}, n::Integer) where {N} = RootOne{N}(r.k * mod(n, N))
Base.inv(r::RootOne{N}) where {N} = RootOne{N}(N - r.k)

Base.literal_pow(::typeof(Base.:^), r::RootOne{N}, ::Val{n}) where {n, N} = RootOne{N}(r.k * mod(n, N))

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

Base.Complex(r::RootOne{2}) = Complex{Int}(r)

function Base.complex(r::RootOne{4})
    k = r.k
    tup =
        if k == 0
            (1, 0)
        elseif k == 1
            (0, 1)
        elseif k == 2
            (-1, 0)
        elseif k == 3
            (0, -1)
        end
    Complex(tup...)
end

"""
    rand([rng=default_rng()], RootOne{N}, [dims...]) where {N}

Return a random element or array of `RootOne{N}`.

# Example
```julia-repl
julia> rand(Omega)
ω₈³

julia> rand(Omega, 3)
3-element Vector{Omega}:
 ω₈⁵
 ω₈²
 ω₈⁵
```
"""
function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{RootOne{N}}) where {N}
    RootOne{N}(rand(rng, 0:N-1))
end

end # module RootOnes
