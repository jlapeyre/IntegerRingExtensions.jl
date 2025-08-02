module RootOnes

import ..Utils: subscript, superscript
import ..Common: sqrt_imaginary, imaginary
import Base: convert, show

export RootOne, RootOne8

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

`N`th roots of unity

`N` should always be a literal or a `const`. Otherwise performance is severely degraded.

The type `const RootOne8 = RootOne{8}` is defined.
`k` will be stored as `mod(k, N)`, which takes values in `(0,...,N-1)`.
# Examples
```jldoctest

julia> RootOne8(1)
RootOne8(1)
```

```jldoctest
julia> RootOne8(-1)
RootOne8(7)

julia> -RootOne8(1) # Unary minus
RootOne8(5)

julia> -RootOne8(5)
RootOne8(1)

julia> -RootOne{9}(1)
ERROR: ArgumentError: InexactError unary minus of type RootOne{9}

julia> RootOne8(1)^8
RootOne8(0)

julia> isone(RootOne8(1)^8)
true

julia> one(RootOne8)
RootOne8(0)

julia> one(RootOne8(5))
RootOne8(0)

julia> RootOne8(2) * RootOne8(3)
RootOne8(5)
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

# This is largely for performance. If N is not constant,
# performance suffers by factor of 100X or more.

"""
    RootOne8(n=1)

Eighth roots of unity
"""
const RootOne8 = RootOne{8}

"""
    imaginary(::Type{RootOne8})

The value of type `RootOne8` that represents the imaginary unit.

This is `RootOne8(2)`.
"""
function imaginary(::Type{RootOne8})
    RootOne8(2)
end

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
    sqrt_imaginary(::Type{RootOne8})

The value of type `RootOne8` that represents the principal square root of the imaginary unit.

This is `RootOne8(1)`.
"""
function sqrt_imaginary(::Type{RootOne8})
    RootOne8(1)
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

function (::Type{T})(r::RootOne{M}) where {M, T <: Real}
    iszero(r.k) && return one(T)
    (iseven(M) && div(M, r.k) == 2) || throw(ArgumentError(lazy"InexactError converting RootOne{$M}($r) to $T"))
    return -one(T)
end

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

Base.float(r::RootOne) = complex(r)
Base.complex(r::RootOne) = Complex(r)
Base.big(r::RootOne) = Complex{BigFloat}(r)
convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N} = Complex{T}(r)

function Complex{T}(r::RootOne{N}) where {T, N}
    kT = convert(T, r.k)
    cispi(2 * kT / N)
end

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

Base.isreal(r::RootOne{N}) where {N} = r.k == 0 || div(N, r.k) == 2
Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k + r2.k)
Base.:/(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k - r2.k)

Base.:^(r::RootOne{N}, n::Integer) where {N} = RootOne{N}(r.k * mod(n, N))
Base.inv(r::RootOne{N}) where {N} = RootOne{N}(N - r.k)

Base.literal_pow(::typeof(Base.:^), r::RootOne{N}, ::Val{n}) where {n, N} = RootOne{N}(r.k * mod(n, N))

function Base.:-(r::RootOne{N}) where {N}
    iseven(N) || throw(ArgumentError(lazy"InexactError unary minus of type RootOne{$N}"))
    RootOne{N}(r.k + div(N, 2))
end

end # module RootOnes
