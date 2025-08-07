module DyadicFractions

import Base: zero, iszero, one, convert, promote_rule, show
import ..Utils: superscript, iszero_strong, isone_strong, greater_than_strong,
    PRETTY, lobit
import ..Common: canonical, mul_half, mul_two, params, conj_root_two
import ..Singletons: InvTwo, InvTwoT, TwoT, Pow
import ILog2

########################
####
#### DyadicFraction
####
########################

"""
    struct DyadicFraction{aT, kT}
        a::aT
        k::kT
    end

Represents the ring `𝔻 = ℤ[½]`.

An instance represents the element `a / 2^k`.

# Examples
```jldoctest
julia> DyadicFraction(5, 3)
5/2³

julia> params(DyadicFraction(5, 3))
(5, 3)

julia> DyadicFraction(5)
5

julia> params(DyadicFraction(5))
(5, 0)
```
"""
struct DyadicFraction{aT<:Integer, kT<:Integer} <: Real
    a::aT
    k::kT
end

"""
    params(d::DyadicFraction)

Return a `Tuple` of the two paramters of `d`: `a` and `k`.

# Example
```jldoctest
julia> x = DyadicFraction(5, 2)
5/2²

julia> params(x)
(5, 2)
```
"""
params(d::DyadicFraction) = (d.a, d.k)

Base.conj(d::DyadicFraction) = d
Base.adjoint(d::DyadicFraction) = d
conj_root_two(d::DyadicFraction) = d
Base.transpose(d::DyadicFraction) = d

Base.:*(::InvTwoT, f::DyadicFraction) = mul_half(f)
Base.:*(::TwoT, f::DyadicFraction) = mul_half(f, -1)

Base.:*(pow::Pow{TwoT}, f::DyadicFraction) = mul_half(f, -pow.n)
Base.:*(pow::Pow{InvTwoT}, f::DyadicFraction) = mul_half(f, pow.n)

"""
    mul_half(f::DyadicFraction{T,V}, n::Integer=1) where {T,V}

Divide `f` by `2^n`.

The value returned is in canonical form.
"""
function mul_half(f::DyadicFraction{T,V}, n::Integer=1) where {T,V}
    mul_two(f, -n)
    # n == 0 && return f
    # if n < 0
    #     canonical(DyadicFraction((T(1) << -n) * f.a, f.k))
    # else
    #     canonical(DyadicFraction(f.a, f.k + V(n)))
    # end
end

function mul_two(f::DyadicFraction{T,V}, n::Integer=1) where {T,V}
    n == 0 && return f
    fr = if n < 0
        DyadicFraction(f.a, f.k + -V(n))
    else
        m = n - f.k
        if m >= 0
            DyadicFraction((T(1) << m) * f.a, zero(f.k))
        else
            DyadicFraction(f.a, -m)
        end
    end
    return canonical(fr)
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

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:BigFloat}
    BigFloat
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:BigInt, V1<:Integer, T<:AbstractFloat}
    BigFloat
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:Base.BitInteger, V1<:Integer, T<:AbstractFloat}
    Float64
end

function promote_rule(::Type{DyadicFraction{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:AbstractFloat}
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

function Base.show(io::IO, ::PRETTY, tup::NTuple{N, <:DyadicFraction}) where {N}
    print(io, "(")
    for (i, d) in enumerate(tup)
        show(io, PRETTY(), d)
        if i < length(tup)
            print(io, ", ")
        end
    end
    print(io, ")")
end

function show(io::IO, ::MIME"text/plain", df::DyadicFraction)
#    print(io, df.a, " / 2^", df.k)
    #    print(io, df.a, " / 2", superscript(df.k))
    if iszero_strong(df.a)
        print(io, zero(df.a))
    else
        if iszero_strong(df.k)
            print(io, df.a)
        else
            print(io, df.a, "/2")
            if !isone_strong(df.k)
               print(io, superscript(df.k))
            end
        end
    end
end

function show(io::IO, ::MIME"text/plain", x::Complex{<:DyadicFraction})
    show(io, PRETTY(), x.re)
    print(io, " + ")
    show(io, PRETTY(), x.im)
end

"""
    canonical(df::DyadicFraction)

Return `df` in canonical form. That is, with `df.k` as small as possible.
"""
function canonical(df::DyadicFraction)
    (num, dexp) = (df.a, df.k)
    iszero(num) && return DyadicFraction(zero(num), zero(dexp))
    n = lobit(num)
    dexp >= n && return DyadicFraction(num >> n, dexp - n)
    DyadicFraction(num >> dexp, zero(dexp))
end

function zero(::Type{DyadicFraction{aT, kT}}) where {aT, kT}
    DyadicFraction(zero(aT), zero(kT))
end

zero(::DyadicFraction{aT, kT}) where {aT, kT} = zero(DyadicFraction{aT, kT})

# Careful. There is more than on way to represent zero.
function iszero(df::DyadicFraction)
    iszero(df.a)
end

function one(::Type{DyadicFraction{aT, kT}}) where {aT, kT}
    DyadicFraction(one(aT), zero(kT))
end

one(::DyadicFraction{aT, kT}) where {aT, kT} = one(DyadicFraction{aT, kT})

# Hmmm, assume it is reduced (canonical) for the moment
function Base.isone(df::DyadicFraction{aT, kT}) where {aT, kT}
    isone_strong(df.a) && iszero_strong(df.k)
end

Base.abs(df::DyadicFraction) = DyadicFraction(abs(df.a), df.k)
Base.abs2(df::DyadicFraction) = df * df
Base.sign(df::DyadicFraction) = sign(df.a)

function convert(::Type{T}, f::DyadicFraction) where {T <: Number}
    if iszero_strong(f.k)
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

function convert(::Type{Rational{T}}, f::DyadicFraction) where {T}
    #  Rational{T}(convert(T, f.a), convert(T, 2)^f.k)
    # shifting is a bit faster than power of two (even though 2 is literal)
    Rational{T}(convert(T, f.a), convert(T, 1) << f.k)
end

DyadicFraction(r::Rational) = convert(DyadicFraction, r)
DyadicFraction(n::Integer) = DyadicFraction(n, zero(n))
DyadicFraction(x::DyadicFraction) = x
DyadicFraction{T,V}(x::DyadicFraction{T,V}) where {T,V} = x

Base.Rational(f::DyadicFraction{aT}) where {aT} = convert(Rational{aT}, f)

function convert(::Type{DyadicFraction}, r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    DyadicFraction(r.num, ILog2.ilog2(r.den))
end

function Complex{Tc}(df::DyadicFraction) where {Tc<:Real}
    T = Complex{Tc}
    T(df.a) * T(1//2)^df.k
end

Base.:(==)(df1::DyadicFraction, df2::DyadicFraction) = _cmp_dyadic(df1, df2, ==)
Base.:(<)(df1::DyadicFraction, df2::DyadicFraction) = _cmp_dyadic(df1, df2, <)
Base.:(<=)(df1::DyadicFraction, df2::DyadicFraction) = _cmp_dyadic(df1, df2, <=)
Base.:(>)(df1::DyadicFraction, df2::DyadicFraction) = _cmp_dyadic(df1, df2, >)
Base.:(>=)(df1::DyadicFraction, df2::DyadicFraction) = _cmp_dyadic(df1, df2, >=)

function _cmp_dyadic(df1::DyadicFraction, df2::DyadicFraction, cmp_func)
    dk = df1.k - df2.k
    if dk > 0
        return cmp_func(df1.a, (1 << dk) * df2.a)
    elseif dk < 0
        return cmp_func((1 << -dk) * df1.a, df2.a)
    else
        return cmp_func(df1.a, df2.a)
    end
end

function Base.:(==)(x::Integer, df::DyadicFraction)
    df1 = canonical(df)
    df1.k == 0 && df1.a == x
end
Base.:(==)(df::DyadicFraction, x::Integer) = x == df


function Base.AbstractFloat(f::DyadicFraction{T,V}) where {T,V}
    W = float(T)
    convert(W, f)
end

#function Base.float(f::DyadicFraction) = convert(Float64, f)

Base.big(f::DyadicFraction) = convert(BigFloat, f)
Base.:*(f1::DyadicFraction, f2::DyadicFraction) = DyadicFraction(f1.a * f2.a, f1.k + f2.k)
Base.:-(f::DyadicFraction) = DyadicFraction(-f.a, f.k)

# We can do one of
# 1. Multiply in numerator and leave unsimplified (canonicalized)
# 2. Mulitply in numerator, then canonicalize
# 3. If n is even, decrement denominator and call again with div(n,2)
# Three seems to be much faster than 2. Even with n = 2^20
function Base.:*(n::Integer, f::DyadicFraction)
    iseven(n) && greater_than_strong(f.k, 0) ?
        (n >> 1) * DyadicFraction(f.a, f.k - 1) :
        DyadicFraction(n * f.a, f.k)
end

Base.:*(f::DyadicFraction, n::Integer) = n * f

Base.:+(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, +)
Base.:-(f1::DyadicFraction, f2::DyadicFraction) = _plus(f1, f2, -)
function _plus(f1::DyadicFraction, f2::DyadicFraction, op)
    (minex, maxex) = minmax(f1.k, f2.k)
    DyadicFraction(op(1 << (f2.k - minex) * f1.a,  1 << (f1.k - minex) * f2.a), maxex)
end

for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(df::DyadicFraction)
        convert($Ti, df)
    end
end

end # module DyadicFractions
