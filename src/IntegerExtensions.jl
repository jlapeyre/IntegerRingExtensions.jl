module IntegerExtensions

import LinearAlgebra

export QuadraticInteger

# Assume no k s.t. D != 4k + 1
# i.e. D - 1 not a multiple of 4
# e.g. D == 5 not allowed
struct QuadraticInteger{D, IntT}
    a::IntT
    b::IntT
end

QuadraticInteger(a::T, b::T, D) where {T <: Integer} = QuadraticInteger{D, T}(a, b)
QuadraticInteger{D}(a::T, b::T) where {T <: Integer, D} = QuadraticInteger{D, T}(a, b)

# Maybe not needed
function Base.copy(q::QuadraticInteger{D}) where D
    QuadraticInteger{D}(q.a, q.b)
end

function Base.convert(::Type{QuadraticInteger{D, IntT}}, q::QuadraticInteger{D}) where {D, IntT}
    QuadraticInteger{D}(IntT(q.a), IntT(q.b))
end

function Base.zero(::Type{QuadraticInteger{D, IntT}}) where {D, IntT}
    QuadraticInteger{D, IntT}(0, 0)
end

Base.zero(q::QuadraticInteger) = zero(typeof(q))

function Base.big(q::QuadraticInteger{D}) where D
    convert(QuadraticInteger{D, BigInt}, q)
end

LinearAlgebra.norm(qi::QuadraticInteger{D}) where D = qi.a^2 - D * qi.b^2

Base.conj(qi::QuadraticInteger{D}) where D = QuadraticInteger{D}(qi.a, -qi.b)

function Base.:*(q1::QuadraticInteger{D}, q2::QuadraticInteger{D}) where D
    QuadraticInteger{D}(q1.a * q2.a + 2 * q1.b * q2.b, q1.a * q2.b + q1.b * q2. a)
end

function Base.:-(q1::QuadraticInteger{D}, q2::QuadraticInteger{D}) where D
    QuadraticInteger{D}(q1.a - q2.a, q1.b - q2.b)
end

function Base.:+(q1::QuadraticInteger{D}, q2::QuadraticInteger{D}) where D
    QuadraticInteger{D}(q1.a + q2.a, q1.b + q2.b)
end

# Is there a better method than this?
function Base.:^(q::QuadraticInteger{D}, n::Integer) where D
    Base.power_by_squaring(q, n)
    # q2 = q
    # for _ in 1:(n - 1)
    #     q2 = q2 * q
    # end
    # return q2
end

end # module IntegerExtensions
