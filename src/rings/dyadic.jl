module Dyadics

import Base: zero, iszero, one, convert, promote_rule, show
import Random
import ..Utils: superscript, iszero_strong, isone_strong, greater_than_strong,
    PRETTY, lobit
import ..Common: canonical, mul_half, mul_two, params, conj_root_two, isunit, invstrict
import ..Singletons: InvTwo, InvTwoT, TwoT, Pow
import ILog2

########################
####
#### Dyadic
####
########################

"""
    struct Dyadic{aT<:Integer, kT<:Integer} <: Real
        a::aT
        k::kT
    end

Represents the ring `𝔻 = ℤ[½]`.

An instance represents the element `a / 2^k`.

# Examples
```jldoctest
julia> Dyadic(5, 3)
5/2³

julia> params(Dyadic(5, 3))
(5, 3)

julia> Dyadic(5)
5

julia> params(Dyadic(5))
(5, 0)
```
"""
struct Dyadic{aT<:Integer, kT<:Integer} <: Real
    function Dyadic{Ta, Tk}(a::Ta, k::Tk) where {Ta, Tk}
        iszero(a) && return new{Ta, Tk}(zero(a), zero(k))
        return new{Ta, Tk}(a, k)
    end
    function Dyadic(a::Ta, k::Tk) where {Ta, Tk}
        iszero(a) && return new{Ta, Tk}(zero(a), zero(k))
        return new{Ta, Tk}(a, k)
    end

    a::aT
    k::kT
end

"""
    params(d::Dyadic)

Return a `Tuple` of the two paramters of `d`: `a` and `k`.

# Example
```jldoctest
julia> x = Dyadic(5, 2)
5/2²

julia> params(x)
(5, 2)
```
"""
params(d::Dyadic) = (d.a, d.k)

Base.conj(d::Dyadic) = d
Base.adjoint(d::Dyadic) = d
conj_root_two(d::Dyadic) = d
Base.transpose(d::Dyadic) = d

Base.:*(::InvTwoT, f::Dyadic) = mul_half(f)
Base.:*(::TwoT, f::Dyadic) = mul_half(f, -1)

Base.:*(pow::Pow{TwoT}, f::Dyadic) = mul_half(f, -pow.n)
Base.:*(pow::Pow{InvTwoT}, f::Dyadic) = mul_half(f, pow.n)

"""
    mul_half(f::Dyadic{T,V}, n::Integer=1) where {T,V}

Divide `f` by `2^n`.

The value returned is in canonical form.
"""
function mul_half(f::Dyadic{T,V}, n::Integer=1) where {T,V}
    mul_two(f, -n)
end

function mul_two(f::Dyadic{T,V}, n::Integer=1) where {T,V}
    n == 0 && return f
    fr = if n < 0
        Dyadic(f.a, f.k + -V(n))
    else
        m = n - f.k
        if m >= 0
            Dyadic((T(1) << m) * f.a, zero(f.k))
        else
            Dyadic(f.a, -m)
        end
    end
    return canonical(fr)
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{Dyadic{T2,V2}}) where {T1,T2,V1,V2}
    T = promote_type(T1, T2)
    V = promote_type(V1, V2)
    Dyadic{T, V}
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:Integer}
    T3 = promote_type(T1, T)
    Dyadic{T3, V1}
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:BigFloat}
    BigFloat
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{T}) where {T1<:BigInt, V1<:Integer, T<:AbstractFloat}
    BigFloat
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{T}) where {T1<:Base.BitInteger, V1<:Integer, T<:AbstractFloat}
    Float64
end

function promote_rule(::Type{Dyadic{T1,V1}}, ::Type{T}) where {T1<:Integer, V1<:Integer, T<:AbstractFloat}
    T3 = promote_type(T1, T)
    Dyadic{T3, V1}
end


function Base.:^(df::Dyadic, n::Integer)
    # Can't copy Int, but really should copy BigInt
    # Need logic for this??
    n == 1 && return df
    return Base.power_by_squaring(df, n)
end

function promote_rule(::Type{Dyadic{Int,Int}}, ::Type{Int})
    Dyadic{Int, Int}
end

function Base.show(io::IO, ::PRETTY, tup::NTuple{N, <:Dyadic}) where {N}
    print(io, "(")
    for (i, d) in enumerate(tup)
        show(io, PRETTY(), d)
        if i < length(tup)
            print(io, ", ")
        end
    end
    print(io, ")")
end

function show(io::IO, ::MIME"text/plain", df::Dyadic)
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

function show(io::IO, ::MIME"text/plain", x::Complex{<:Dyadic})
    show(io, PRETTY(), x.re)
    print(io, " + ")
    show(io, PRETTY(), x.im)
end

"""
    canonical(df::Dyadic)

Return `df` in canonical form. That is, with `df.k` as small as possible.
"""
function canonical(df::Dyadic)
    (num, dexp) = (df.a, df.k)
    iszero(num) && return Dyadic(zero(num), zero(dexp))
    n = lobit(num)
    dexp >= n && return Dyadic(num >> n, dexp - n)
    Dyadic(num >> dexp, zero(dexp))
end

function zero(::Type{Dyadic{aT, kT}}) where {aT, kT}
    Dyadic(zero(aT), zero(kT))
end

zero(::Dyadic{aT, kT}) where {aT, kT} = zero(Dyadic{aT, kT})

# Careful. There is more than on way to represent zero.
function iszero(df::Dyadic)
    iszero(df.a)
end

function one(::Type{Dyadic{aT, kT}}) where {aT, kT}
    Dyadic(one(aT), zero(kT))
end

one(::Dyadic{aT, kT}) where {aT, kT} = one(Dyadic{aT, kT})

function Base.isone(df::Dyadic{aT, kT}) where {aT, kT}
    df1 = canonical(df)
    isone_strong(df1.a) && iszero_strong(df1.k)
end

function Base.isinteger(df::Dyadic)
    return df.k <= 0
end

function isunit(df::Dyadic)
    return ispow2(df.a)
end

function invstrict(df::Dyadic)
    isunit(df) || throw(ArgumentError(lazy"$df has no inverse of type $(typeof(df))"))
    return typeof(df)(1, ILog2.ilog2(df.a) - df.k)
end

Base.abs(df::Dyadic) = Dyadic(abs(df.a), df.k)
Base.abs2(df::Dyadic) = df * df
Base.sign(df::Dyadic) = sign(df.a)

## According to the Julia manual, we are doing this backward.
## We should define the constructors. And then the convert methods when
## neccessary. Becuase `convert` is rather magic.

function convert(::Type{T}, f::Dyadic) where {T <: Number}
    if iszero_strong(f.k)
        return convert(T, f.a)
    end
    (a, two) = promote(f.a, 2)
    convert(T, a) / convert(T, two)^f.k
end

function convert(::Type{Float64}, f::Dyadic)
    if iszero(f.k)
        return convert(Float64, f.a)
    end
    convert(Float64, f.a) * 0.5^f.k
end

function convert(::Type{Dyadic{T,V}}, f::Dyadic) where {T, V}
    Dyadic{T, V}(T(f.a), V(f.k))
end

function convert(::Type{T}, f::Dyadic) where {T <: Integer}
    iszero(f.k) && return convert(T, f.a)
    cf = canonical(f)
    cf.k == 0 || throw(ArgumentError(lazy"Inexact error converting $f to $T"))
    convert(T, cf.a)
end

function convert(::Type{Dyadic}, n::Integer)
    Dyadic(n, zero(n))
end

convert(::Type{Dyadic{T, K}}, n::T1) where {T <: Integer, T1 <: Integer, K} = Dyadic{T, K}(n)

function Dyadic{T, K}(n::T1)  where {T <: Integer, T1 <: Integer, K}
    Dyadic{T, K}(convert(T, n), 0)
end

function convert(::Type{Rational{T}}, f::Dyadic) where {T}
    #  Rational{T}(convert(T, f.a), convert(T, 2)^f.k)
    # shifting is a bit faster than power of two (even though 2 is literal)
    Rational{T}(convert(T, f.a), convert(T, 1) << f.k)
end

function Dyadic(r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    Dyadic(r.num, ILog2.ilog2(r.den))
end

function Dyadic{T, V}(r::Rational) where {T<:Integer, V<:Integer}
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    Dyadic{T,V}(r.num, ILog2.ilog2(r.den))
end

#Dyadic(r::Rational) = convert(Dyadic, r)
Dyadic(n::Integer) = Dyadic(n, zero(n))
Dyadic(x::Dyadic) = x
Dyadic{T,V}(x::Dyadic{T,V}) where {T,V} = x

Base.Rational(f::Dyadic{aT}) where {aT} = convert(Rational{aT}, f)

function convert(::Type{Dyadic}, r::Rational)
    ispow2(r.den) || throw(ArgumentError(lazy"denominator $(r.den) not a power of 2"))
    Dyadic(r.num, ILog2.ilog2(r.den))
end

function Complex{Tc}(df::Dyadic) where {Tc<:Real}
    T = Complex{Tc}
    T(df.a) * T(1//2)^df.k
end

Base.:(==)(df1::Dyadic, df2::Dyadic) = _cmp_dyadic(df1, df2, ==)
Base.:(<)(df1::Dyadic, df2::Dyadic) = _cmp_dyadic(df1, df2, <)
Base.:(<=)(df1::Dyadic, df2::Dyadic) = _cmp_dyadic(df1, df2, <=)
Base.:(>)(df1::Dyadic, df2::Dyadic) = _cmp_dyadic(df1, df2, >)
Base.:(>=)(df1::Dyadic, df2::Dyadic) = _cmp_dyadic(df1, df2, >=)

function _cmp_dyadic(df1::Dyadic, df2::Dyadic, cmp_func)
    dk = df1.k - df2.k
    if dk > 0
        return cmp_func(df1.a, (1 << dk) * df2.a)
    elseif dk < 0
        return cmp_func((1 << -dk) * df1.a, df2.a)
    else
        return cmp_func(df1.a, df2.a)
    end
end

function Base.:(==)(x::Integer, df::Dyadic)
    df1 = canonical(df)
    df1.k == 0 && df1.a == x
end
Base.:(==)(df::Dyadic, x::Integer) = x == df


function Base.AbstractFloat(f::Dyadic{T,V}) where {T,V}
    W = float(T)
    convert(W, f)
end

#function Base.float(f::Dyadic) = convert(Float64, f)

Base.big(f::Dyadic) = convert(BigFloat, f)
Base.:*(f1::Dyadic, f2::Dyadic) = Dyadic(f1.a * f2.a, f1.k + f2.k)
Base.:-(f::Dyadic) = Dyadic(-f.a, f.k)

# We can do one of
# 1. Multiply in numerator and leave unsimplified (canonicalized)
# 2. Mulitply in numerator, then canonicalize
# 3. If n is even, decrement denominator and call again with div(n,2)
# Three seems to be much faster than 2. Even with n = 2^20
function Base.:*(n::Integer, f::Dyadic)
    iseven(n) && greater_than_strong(f.k, 0) ?
        (n >> 1) * Dyadic(f.a, f.k - 1) :
        Dyadic(n * f.a, f.k)
end

Base.:*(f::Dyadic, n::Integer) = n * f

Base.:+(f1::Dyadic, f2::Dyadic) = _plus(f1, f2, +)
Base.:-(f1::Dyadic, f2::Dyadic) = _plus(f1, f2, -)
function _plus(f1::Dyadic, f2::Dyadic, op)
    (minex, maxex) = minmax(f1.k, f2.k)
    Dyadic(op(1 << (f2.k - minex) * f1.a,  1 << (f1.k - minex) * f2.a), maxex)
end

"""
    struct DyadicSample{T,V}

# Example

```julia-repl
julia> rand(DyadicSample(1:10^3, 1:10))
509/2⁷

julia> rand(DyadicSample(-3:3, 1:5), 3)
3-element Vector{Dyadic{Int64, Int64}}:
  2/2²
 -3/2⁵
 -3/2
```
"""
struct DyadicSample{T,V}
    asamp::T
    ksamp::V
end

function Base.eltype(::Type{DyadicSample{T,V}}) where {T <: AbstractArray, V <: AbstractArray}
    tt = eltype(T)
    vv = eltype(V)
    Dyadic{tt, vv}
end

function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerTrivial{<:DyadicSample})
    Dyadic(rand(rng, s[].asamp), rand(rng, s[].ksamp))
end

# This is really backward.
# We should define the constructor functions.
# Then, when neccessary, define the convert functions.
for Ti in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
    @eval function (::Type{Base.$Ti})(df::Dyadic)
        convert($Ti, df)
    end
end

end # module Dyadics
