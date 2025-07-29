module Matrices

struct Matrix2x2{T} <: AbstractMatrix{T}
    data::NTuple{4, T}
end

# Hmmm. Depends on whether type T is copyable
# In any case, a Tuple is not.
# Base.copy(m::Matrix2x2) = Matrix2x2(copy(m.data))

function Matrix2x2(a, b, c, d)
    Matrix2x2((a, b, c, d))
end

Base.size(::Matrix2x2) = (2, 2)

Base.IndexStyle(::Type{<:Matrix2x2}) = IndexLinear()

function Base.one(::Type{Matrix2x2{T}}) where T
    Matrix2x2(one(T), zero(T), zero(T), one(T))
end

Base.one(::Matrix2x2{T}) where {T} = one(Matrix2x2{T})

@inline function Base.getindex(m::Matrix2x2, i::Integer)
    @boundscheck checkbounds(m, i)
    return @inbounds m.data[i]
end

@inline function _mul2(m1, m2)
    (a1, b1, c1, d1) = m1.data
    (a2, b2, c2, d2) = m2.data
    Matrix2x2(a1*a2 + c1*b2, a1*c2+c1*d2, b1*a2+d1*b2, b1*c2 + d1*d2)
end

@inline function _add2(m1, m2, op=+)
    (a1, b1, c1, d1) = m1.data
    (a2, b2, c2, d2) = m2.data
    Matrix2x2(op(a1, a2), op(b1, b2), op(c1, c2), op(d1,d2))
end

Base.:*(m1::Matrix2x2, m2::Matrix2x2) = _mul2(m1, m2)
Base.:+(m1::Matrix2x2, m2::Matrix2x2) = _add2(m1, m2)
Base.:-(m1::Matrix2x2, m2::Matrix2x2) = _add2(m1, m2, -)

function power(m::Matrix2x2, n::Integer)
    n == 0 && return one(m)
    n == 1 && return m  # If elements of m are mutable, this is problematic
    n == 2 && return m * m
    return Base.power_by_squaring(m, n)
end

# This is only called if `n` is not literal at the call site.
function Base.:^(m::Matrix2x2, n::Integer)
    power(m, n)
end

Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{0}) = one(m)
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{1}) = m # If elements of m are mutable, this is problematic
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{2}) = m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{3}) = m * m * m
Base.literal_pow(::typeof(Base.:^), m::Matrix2x2, ::Val{4}) = (m * m) * (m * m)

function Matrix2x2{T}(m::Matrix2x2) where {T}
    @inbounds Matrix2x2(convert(T, m[1]), convert(T, m[2]), convert(T, m[3]), convert(T, m[4]))
end

Base.convert(::Type{Matrix2x2{T}}, m::Matrix2x2) where {T} = Matrix2x2{T}(m)

end # module Matrices
