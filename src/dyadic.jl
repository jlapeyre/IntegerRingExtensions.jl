module DyadicFractions

import Base: zero, iszero, one, convert, promote_rule, show
import ..Utils: superscript
import ..Common: canonical

########################
####
#### DyadicFraction
####
########################

"""
    DyadicFraction{aT, kT}

Represents the ring `𝔻 = ℤ[½]`.

`aT` is the type of the numerator.
`bT` is the type of the exponent on `2` in the denominator.

# Examples
```jldoctest
julia> DyadicFraction(5, 3)
5/2³
```
"""
struct DyadicFraction{aT, kT}
    a::aT
    k::kT
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{DyadicFraction{T2,V2}}) where {T1,T2,V1,V2}
    T = promote_type(T1, T2)
    V = promote_type(V1, V2)
    DyadicFraction{T, V}
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:Integer}
    T3 = promote_type(T1, T)
    DyadicFraction{T3, V1}
end

function Base.:^(df::DyadicFraction, n::Integer)
    # Can't copy Int, but really should copy BigInt
    # Need logic for this??
    n == 1 && return df
    return Base.power_by_squaring(df, n)
end

function promote_rule(::Type{DyadicFraction{Int,Int}}, ::Type{Int})
    DyadicFraction{Int, Int}
end

function show(io::IO, ::MIME"text/plain", df::DyadicFraction)
#    print(io, df.a, " / 2^", df.k)
    #    print(io, df.a, " / 2", superscript(df.k))
    if iszero(df.a)
        print(io, zero(df.a))
    else
        if iszero(df.k)
            print(io, df.a)
        else
            print(io, df.a, "/2")
            if !isone(df.k)
               print(io, superscript(df.k))
            end
        end
    end
end

"""
    canonical(df::DyadicFraction)

Return `df` in canonical form. That is, with `df.k` as small as possible.
"""
function canonical(df::DyadicFraction)
    iszero(df.a) && return DyadicFraction(df.a, zero(df.k))
    (num, dexp) = (df.a, df.k)
    while true
        (isodd(num) || dexp < 1) && return DyadicFraction(num, dexp)
        num = div(num, 2)
        dexp = dexp - 1
    end
end

# Routine above is faster and simpler than this one
# function oldcanonical(df::DyadicFraction)
#     iszero(df.a) && return DyadicFraction(df.a, zero(df.k))
#     (c, m) = factortwos(df.a)
#     iszero(c) && return df
#     k = df.k
#     if c > k
#         return typeof(df)(2^(c-k)*m, 0)
#     elseif c < k
#         return typeof(df)(m, k - c)
#     else
#         return typeof(df)(m, 0)
#     end
# end

# Unused now
# """
#     factortwos(n)
#
# Factor `n` as `m * 2^c` and return `(c, m)`.
# """
# function factortwos(n)
#     n == 0 && return (0, 0)
#     c = 0
#     m = n
#     while true
#         isodd(m) && return (c, m)
#         m = div(m, 2)
#         c += 1
#     end
# end

function zero(::Type{DyadicFraction{aT, bT}}) where {aT, bT}
    DyadicFraction(zero(aT), zero(bT))
end

zero(::DyadicFraction{aT, bT}) where {aT, bT} = zero(DyadicFraction{aT, bT})

# Careful. There is more than on way to represent zero.
function iszero(df::DyadicFraction)
    iszero(df.a)
end

function one(::Type{DyadicFraction{aT, bT}}) where {aT, bT}
    DyadicFraction(one(aT), zero(bT))
end

one(::DyadicFraction{aT, bT}) where {aT, bT} = one(DyadicFraction{aT, bT})

function convert(::Type{T}, f::DyadicFraction) where {T}
    if iszero(f.k)
        return convert(T, f.a)
    end
    convert(T, f.a) / convert(T, 2)^f.k
end

function convert(::Type{Float64}, f::DyadicFraction)
    if iszero(f.k)
        return convert(Float64, f.a)
    end
    convert(Float64, f.a) * 0.5^f.k
end

function convert(::Type{DyadicFraction{T,V}}, f::DyadicFraction) where {T, V}
    DyadicFraction{T, V}(T(f.a), V(f.k))
end

function convert(::Type{T}, f::DyadicFraction) where {T <: Integer}
    iszero(f.k) && return convert(T, f.a)
    cf = canonical(f)
    cf.k == 0 || throw(ArgumentError(lazy"Inexact error converting $f to $T"))
    convert(T, cf.a)
end

function convert(::Type{DyadicFraction}, n::Integer)
    DyadicFraction(n, zero(n))
end

convert(::Type{DyadicFraction{T, K}}, n::T1) where {T <: Integer, T1 <: Integer, K} = DyadicFraction{T, K}(n)

function DyadicFraction{T, K}(n::T1)  where {T <: Integer, T1 <: Integer, K}
    DyadicFraction{T, K}(convert(T, n), 0)
end

convert(::Type{Rational{T}}, f::DyadicFraction) where {T} =
    Rational{T}(convert(T, f.a), convert(T, 2)^f.k)

DyadicFraction(r::Rational) = convert(DyadicFraction, r)
DyadicFraction(n::Integer) = DyadicFraction(n, zero(n))
DyadicFraction(x::DyadicFraction) = x
DyadicFraction{T,V}(x::DyadicFraction{T,V}) where {T,V} = x

Base.Rational(f::DyadicFraction{aT}) where {aT} = convert(Rational{aT}, f)

function convert(::Type{DyadicFraction}, r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    DyadicFraction(r.num, ilog2(r.den))
end

function Complex{Tc}(df::DyadicFraction) where {Tc}
    T = Complex{Tc}
    T(df.a) * T(1//2)^df.k
end

function Base.:(==)(df1::DyadicFraction, df2::DyadicFraction)
    if df1.k > df2.k
        return df1.a == 2^(df1.k - df2.k) * df2.a
    elseif df1.k < df2.k
        return df2.a == 2^(df2.k - df1.k) * df1.a
    else
        return df1.a == df2.a
    end
end


Base.float(f::DyadicFraction) = convert(Float64, f)
Base.big(f::DyadicFraction) = convert(BigFloat, f)
Base.:*(f1::DyadicFraction, f2::DyadicFraction) = DyadicFraction(f1.a * f2.a, f1.k + f2.k)
Base.:-(f::DyadicFraction) = DyadicFraction(-f.a, f.k)

# We can do one of
# 1. Multiply in numerator and leave unsimplified (canonicalized)
# 2. Mulitply in numerator, then canonicalize
# 3. If n is even, decrement denominator and call again with div(n,2)
# Three seems to be much faster than 2. Even with n = 2^20
function Base.:*(n::Integer, f::DyadicFraction)
    iseven(n) && f.k > 0 ?
        div(n, 2) * DyadicFraction(f.a, f.k - 1) :
        DyadicFraction(n * f.a, f.k)
end

Base.:*(f::DyadicFraction, n::Integer) = n * f

Base.:+(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, +)
Base.:-(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, -)
function _plus(f1::DyadicFraction, f2::DyadicFraction, op)
    (minex, maxex) = minmax(f1.k, f2.k)
    DyadicFraction(op(2^(f2.k - minex) * f1.a,  2^(f1.k - minex) * f2.a), maxex)
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(df::DyadicFraction)
        convert($Ti, df)
    end
end

end # module DyadicFractions
