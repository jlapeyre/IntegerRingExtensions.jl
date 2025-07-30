module RootOnes

import ..IntegerExtensions.Utils: subscript, superscript
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
    sqrt_imaginary(::Type{RootOne8})

The value of type `RootOne8` that represents the principal square root of the imaginary unit.

This is `RootOne8(1)`.
"""
function sqrt_imaginary(::Type{RootOne8})
    RootOne8(1)
end

function show(io::IO, ::MIME"text/plain", r::RootOne{N}) where {N}
    print(io, "ω" * subscript(N))
    isone(r.k) || print(io, superscript(r.k))
end

convert(::Type{Complex{T}}, r::RootOne{N}) where {T, N} = Complex{T}(r)

function Complex{T}(r::RootOne{N}) where {T, N}
    kT = convert(T, r.k)
    cispi(2 * kT / N)
end

"""
    RootOne(k, N)

Construct `RootOne{N}(k)`.
"""
RootOne(k, N) = RootOne{N}(k)
Base.one(::Type{RootOne{N}}) where {N} = RootOne{N}(0)
Base.one(x::RootOne) = one(typeof(x))

Base.float(r::RootOne) = complex(r)
Base.complex(r::RootOne) = convert(ComplexF64, r)
Base.big(r::RootOne) = convert(Complex{BigFloat}, r)
# TODO: angle should agree with floating point angle: in [-pi, pi] rather than [0, 2pi]

function Base.angle(r::RootOne{N}) where {N}
    k = r.k <= N/2 ? r.k : r.k - N
    2 * pi * k / N
end

function Base.angle(::Type{T}, r::RootOne{N}) where {N, T}
    k = r.k <= N/2 ? r.k : r.k - N
    2 * T(pi) * T(k) / T(N)
end

Base.conj(r::RootOne) = inv(r)
Base.abs2(r::RootOne) = 1
Base.abs(r::RootOne) = 1
Base.real(r::RootOne{N}) where {N} = cospi(2 * r.k / N)
Base.imag(r::RootOne{N}) where {N} = sinpi(2 * r.k / N)
Base.isreal(r::RootOne{N}) where {N} = r.k == 0 || div(N, r.k) == 2

Base.:*(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k + r2.k)
Base.:/(r1::RootOne{N}, r2::RootOne{N}) where {N} = RootOne{N}(r1.k - r2.k)

Base.:^(r::RootOne{N}, n::Integer) where {N} = RootOne{N}(r.k * n)
Base.inv(r::RootOne{N}) where {N} = RootOne{N}(N - r.k)

function Base.:-(r::RootOne{N}) where {N}
    iseven(N) || throw(ArgumentError(lazy"InexactError unary minus of type RootOne{$N}"))
    RootOne{N}(r.k + div(N, 2))
end

end # module RootOnes
