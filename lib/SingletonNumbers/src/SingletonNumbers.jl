"""
    module SingletonNumbers

This module implements several types, each representing a single number.

They types are intended for use in arithmetic expressions involving elements of ring-extensions of the
integers. The types are meant to be used in  expressions that are evaluated at compile time. Otherwise, dynamic dispatch will seriously
degrade performance.

### Singleton number types
The singleton types are subtypes of the abstract type [`SingleNum`](@ref).

- [`Zero`](@ref)
- [`One`](@ref)
- [`Two`](@ref)
- [`InvTwo`](@ref) Reciprocal of two.
- [`RootTwo`](@ref) Square root of two.
- [`InvRootTwo`](@ref) Reciprocal of square root of two.
- [`Imag`](@ref) Imaginary unit, 𝕚.
- [`RootImag`](@ref) Square root of imaginary unit, √𝕚.

### `struct Pow`
- [`Pow`](@ref) -- `Pow{T}(n)` represents a singleton number of type `T` raised to the power `n`.
"""
module SingletonNumbers

import Base: show, inv, sqrt, isone, iszero, isinteger, iseven, isreal
import RingExtensionsUtils: PRETTY, superscript

export RootTwo, InvRootTwo, Imag, 𝕚,  RootImag, Two, 𝟚, InvTwo, 𝟚⁻¹,
    One, Zero, 𝟙, 𝟘

export ZeroT, OneT, TwoT, InvTwoT, RootTwoT, InvRootTwoT, ImagT, RootImagT

import RingExtensionsCommon: isrational
export isrational

export Pow

###
### Zero
###

# Should this be a Number? What is a `Number`?
"""
    abstract type SingleNum

An abstract type for reprenting single numbers as singleton types.
Each subtype implements a method for `Base` functions
`isinteger`, `iszero`, `isone`, `iseven`, `isreal`, `inv`.

Each type implements `RingExtensionsCommon.isrational`.
"""
abstract type SingleNum end

show(io::IO, x::SingleNum) = show(io, PRETTY(), x)

const NumSingleNum = Union{Number, SingleNum}

struct ZeroT <: SingleNum
end
const Zero = ZeroT()
const 𝟘 = Zero
show(io::IO, ::PRETTY, ::ZeroT) = print(io, "𝟘")

iszero(::ZeroT) = true
isone(::ZeroT) = false
iseven(::ZeroT) = true
isinteger(::ZeroT) = true
isrational(::ZeroT) = true
isreal(::ZeroT) = true
inv(::ZeroT) = Inf

"""
    Zero
    𝟘

Represents the additive identity, 𝟘, in the ring of integers.

`Zero` also represents the additive identity in any ring or field that extends the integers.

### Examples
```jldoctest
julia> Zero
𝟘

julia> 𝟘
𝟘

julia> TwoT() === Two
true

julia> Tuple(f(Zero) for f in (iszero, isone, iseven, isinteger, isrational, isreal, inv))
(true, false, true, true, true, true, Inf)

julia> Tuple(Zero * x for x in (Zero, One, Two, InvTwo, RootTwo, RootTwo, InvRootTwo, Imag, RootImag))
(𝟘, 𝟘, 𝟘, 𝟘, 𝟘, 𝟘, 𝟘, 𝟘, 𝟘)

julia> Tuple(Zero + x for x in (Zero, One, Two, InvTwo, RootTwo, RootTwo, InvRootTwo, Imag, RootImag))
(𝟘, 𝟙, 𝟚, 𝟚⁻¹, √𝟚, √𝟚, √𝟚⁻¹, 𝕚, √𝕚)
```
"""
Zero

@doc (@doc Zero) 𝟘

###
### One
###

struct OneT <: SingleNum
end
const One = OneT()
const 𝟙 = One
show(io::IO, ::PRETTY, ::OneT) = print(io, "𝟙")

isrational(::OneT) = true
isinteger(::OneT) = true
iszero(::OneT) = false
isone(::OneT) = true
iseven(::OneT) = false
isreal(::OneT) = true

Base.zero(n::SingleNum) = Zero
Base.one(n::SingleNum) = One

"""
    One
    𝟙

Represents the multiplicative identity, 𝟙, in the ring of integers.
"""
One

@doc (@doc One) 𝟙

###
### Two
###

struct TwoT  <: SingleNum
end
const Two = TwoT()
const 𝟚 = Two
show(io::IO, ::PRETTY, ::TwoT) = print(io, "𝟚")

isrational(::TwoT) = true
isinteger(::TwoT) = true
iszero(::TwoT) = false
isone(::TwoT) = false
iseven(::TwoT) = true
isreal(::TwoT) = true

"""
    Two
    𝟚

Represents the integer two: 𝟚
"""
Two

@doc (@doc Two) 𝟚

###
### InvTwo
###

struct InvTwoT <: SingleNum
end
const InvTwo = InvTwoT()
const 𝟚⁻¹ = InvTwo
show(io::IO, ::PRETTY, ::InvTwoT) = print(io, "𝟚⁻¹")
inv(::TwoT) = InvTwo
inv(::InvTwoT) = Two
isrational(::InvTwoT) = true
isinteger(::InvTwoT) = false
iszero(::InvTwoT) = false
isone(::InvTwoT) = false
iseven(::InvTwoT) = false
isreal(::InvTwoT) = true

"""
    InvTwo
    𝟚⁻¹

Represents the multiplicative inverse of the integer two: 𝟚⁻¹
"""
InvTwo

@doc (@doc InvTwo) 𝟚⁻¹

###
### RootTwo
###

struct RootTwoT <: SingleNum
end
const RootTwo = RootTwoT()

show(io::IO, ::PRETTY, ::RootTwoT) = print(io, "√𝟚")
sqrt(::TwoT) = RootTwo
isrational(::RootTwoT) = false
isinteger(::RootTwoT) = false
iszero(::RootTwoT) = false
isone(::RootTwoT) = false
iseven(::RootTwoT) = false
isreal(::RootTwoT) = true

"""
    RootTwo

Represents the square root of two: √𝟚
"""
RootTwo

###
### InvRootTwo
###

struct InvRootTwoT <: SingleNum
end
const InvRootTwo = InvRootTwoT()
show(io::IO, ::PRETTY, ::InvRootTwoT) = print(io, "√𝟚⁻¹")
sqrt(::InvTwoT) = InvRootTwo
inv(::InvRootTwoT) = RootTwo
inv(::RootTwoT) = InvRootTwo
isrational(::InvRootTwoT) = false
isinteger(::InvRootTwoT) = false
iszero(::InvRootTwoT) = false
isone(::InvRootTwoT) = false
iseven(::InvRootTwoT) = false
isreal(::InvRootTwoT) = true

"""
    InvRootTwo

Represents the reciprocoal of the square root of two: √𝟚⁻¹
"""
InvRootTwo

###
### Imag
###

struct ImagT <: SingleNum
end
const Imag = ImagT()
const 𝕚 = Imag
show(io::IO, ::PRETTY, ::ImagT) = print(io, "𝕚")

# Not in Julia Base or stdlib. I define this to include Gaussian rationals
# This should be evaluated based on utility.
isrational(::ImagT) = true

# Julia says imaginary numbers can't be integers. So no Gaussian integers
isinteger(::ImagT) = false
iszero(::ImagT) = false
isone(::ImagT) = false
iseven(::ImagT) = false
isreal(::ImagT) = false

"""
    Imag
    𝕚

Represents the imaginary unit: 𝕚

### Examples
```jldoctest
julia> 𝕚 * 𝕚
-1

julia> 𝕚 + 𝕚
ERROR: MethodError: no method matching +(::ImagT, ::ImagT)
```
"""
Imag

@doc (@doc Imag) 𝕚

###
### RootImag
###

struct RootImagT <: SingleNum
end

const RootImag = RootImagT()
show(io::IO, ::PRETTY, ::RootImagT) = print(io, "√𝕚")
sqrt(::ImagT) = RootImag
isrational(::RootImagT) = false
isinteger(::RootImagT) = false
iszero(::RootImagT) = false
isone(::RootImagT) = false
iseven(::RootImagT) = false
isreal(::RootImagT) = false

"""
    RootImag

Represents the principal square root of the imaginary unit: √𝕚
This is also the principal eight root of one.
"""
RootImag

function Base.conj(x::SingleNum)
    isreal(x) || throw(ArgumentError(lazy"Unsupported"))
    x
end

Base.adjoint(x::SingleNum) = conj(x)

"""
    struct Pow{T<:SingleNum} <: Number
        n::Int
    end


Represents an instance of type `T` raised to the power `n`.
"""
struct Pow{T<:SingleNum} <: Number
    n::Int

    function Pow{T}(n::Int) where T<:SingleNum
        new{T}(n)
    end
end

function show(io::IO, ::PRETTY, p::Pow{T}) where {T}
    show(io, PRETTY(), T())
    print(io, superscript(p.n))
end

function Base.conj(x::Pow{T}) where {T <: SingleNum}
    isreal(T()) || throw(ArgumentError(lazy"Unsupported"))
    x
end

Base.adjoint(x::Pow) = conj(x)

##
## These need testing
##

###
### Write down, somewhere in this monorepo, our understanding of Julia's SOP for
### conversion, etc. We need a brief synopsis to refer to.
### This includes the difference between `float` and `AbstractFloat`.
### Eg. `AbstractFloat` converts to the "best" type `T` s.t. `T <: AbstractFloat`
### OTOH, `float(x)` converts `x` to a type floating-point type, more broadly construed.
### For example `ComplexF64`. It finds the best (narrowest in some sense maybe) type that
### can represent `x` as a broader floating point type.
### So: if `x` can be represented exactly-ish as `T <: AbstractFloat`, then
### Implement just `AbstractFloat`. And use the fallback method for `float` to call this.
### Else, write a special method for `float`.
###
### Some issues. Julia has no abstract Complex type (everyone wishes there were)

function Base.AbstractFloat(p::Pow{T}) where {T}
    AbstractFloat(T())^p.n
end

function Base.Float64(p::Pow{T}) where {T}
    AbstractFloat(T())^p.n
end

# TODO: finish Complex, complex, etc.
function Base.Complex{V}(p::Pow{T}) where {T, V<:Real}
    Complex{V}(T())^p.n
end

function Base.Complex(p::Pow{T}) where {T}
    Complex(T())^p.n
end

function Base.Rational(p::Pow{T}) where {T}
    Rational(T())^p.n
end

function Base.Rational{V}(p::Pow{T}) where {V, T}
    Rational{V}(T())^p.n
end

function show(io::IO, ::PRETTY, p::Pow{InvRootTwoT})
    show(io, PRETTY(), RootTwo)
    print(io, superscript(-p.n))
end

Base.length(::Pow{<:SingleNum}) = 1
Base.iterate(p::Pow{<:SingleNum}) = (p, nothing)
Base.iterate(p::Pow{<:SingleNum}, ::Any) = nothing

Base.length(::SingleNum) = 1
Base.iterate(p::SingleNum) = (p, nothing)
Base.iterate(p::SingleNum, ::Any) = nothing

Base.size(::SingleNum) = ()
Base.size(::Pow{<:SingleNum}) = ()

Base.isinteger(p::Pow{TwoT}) = p.n >= 0
Base.isinteger(p::Pow{InvTwoT}) = p.n <= 0
Base.isinteger(p::Pow{RootTwoT}) = p.n >= 0 && iseven(p.n)
Base.isinteger(p::Pow{InvRootTwoT}) = p.n <= 0 && iseven(p.n)
Base.isinteger(p::Pow{ZeroT}) = p.n >= 0

isrational(p::Pow{TwoT}) = true
isrational(p::Pow{InvTwoT}) = true
isrational(p::Pow{RootTwoT}) = iseven(p.n)
isrational(p::Pow{InvRootTwoT}) = iseven(p.n)
isrational(p::Pow{ZeroT}) = p.n >= 0

###
### Conversion and arithmetic
###

"""
    canconvert(::Type{T}, ::Type{V})::Bool
    canconvert(obj::T, ::Type{V})::Bool

Returns `true` if we can convert `obj::T` with the call `V(obj)`.

The result should be an object of type `<: V`.
"""
function canconvert end

Base.convert(::Type{T}, obj::SingleNum) where {T<:Number} = T(obj)

# If called on an object, call again, on type of object
canconvert(obj::SingleNum, ::Type{V}) where {V} = canconvert(typeof(obj), V)

###
### <: SingleNum
###

## TODO: Re-evaluate whether the methods in this section are necessary, or useful.
## I had to add particular cases due to dispatch ambiguities for cases when these more general methods were intended
## to apply.

function canconvert(::Type{ST}, ::Type{T}) where {ST<:SingleNum, T <: AbstractFloat}
    canconvert(ST, Rational{Int})
end

# Almost no type can be converted to `Bool`. Maybe special case `One` and `Zero`.
# function canconvert(::Type{ST}, ::Type{T}) where {ST<:SingleNum, T<:Bool}
#     false
# end

function canconvert(::Type{ST}, ::Type{T}) where {ST <: Union{OneT, TwoT}, T<:Bool}
    true
end

# If I can convert to `T<:Real`, then I can convert to `Complex{T}`
# i.e. If `Real(obj)` works, then `Complex(obj)` works.
# If the method is missing, this is a violation of the (informal) interface.
function canconvert(::Type{ST}, ::Type{Complex{T}}) where {ST<:SingleNum, T <: Real}
    canconvert(ST, T)
end

# If I can convert to `T<:Integer`, then I can convert to `Rational{T}`
# If `Integer(obj)` works, then `Rational(obj)` must.
function canconvert(::Type{ST}, ::Type{Rational{T}}) where {ST<:SingleNum, T<:Integer}
    canconvert(ST, T)
end

# If I can convert with `Integer(obj)`, then I can convert with `Rational(obj)`.
function canconvert(::Type{ST}, ::Type{Rational}) where {ST<:SingleNum}
    canconvert(ST, Integer)
end

# If `AbstractFloat(obj)` works, then `Complex(obj)` must.
function canconvert(::Type{ST}, ::Type{Complex}) where {ST<:SingleNum}
    canconvert(ST, AbstractFloat)
end

###
### One
###

canconvert(::Type{OneT}, ::Type{T}) where {T <: Integer} = true
canconvert(::Type{OneT}, ::Type{Real}) = true
# Method `Bool(One)` is written explicitly, not via `canconvert`.
canconvert(::Type{OneT}, ::Type{T}) where T<:Bool = true

###
### Two
###

canconvert(::Type{TwoT}, ::Type{T}) where {T <: Integer} = true
canconvert(::Type{TwoT}, ::Type{T}) where T<:Bool = false
canconvert(::Type{TwoT}, ::Type{Real}) = true

for ST in (:RootTwoT, :InvRootTwoT)
    @eval canconvert(::Type{$ST}, ::Type{T}) where {T <: Integer} = false
    @eval canconvert(::Type{$ST}, ::Type{T}) where {T <: AbstractFloat} = true
    @eval canconvert(::Type{$ST}, ::Type{Real}) = true
end

###
### InvTwo
###

canconvert(::Type{InvTwoT}, ::Type{T}) where {T <: Integer} = false
canconvert(::Type{InvTwoT}, ::Type{Rational}) = true
canconvert(::Type{InvTwoT}, ::Type{Real}) = true

function canconvert(::Type{InvTwoT}, ::Type{Rational{T}}) where {T<:Integer}
    canconvert(Two, T)
end

canconvert(::Type{InvTwoT}, ::Type{T}) where {T <: Rational{Bool}} = false

###
### Imag
###

canconvert(::Type{ImagT}, ::Type{T}) where {T <: Real} = false
canconvert(::Type{ImagT}, ::Type{T}) where {T <: Integer} = false
canconvert(::Type{ImagT}, ::Type{Rational}) = false
canconvert(::Type{ImagT}, ::Type{Rational{T}}) where {T <: Integer} = false
canconvert(::Type{ImagT}, ::Type{T}) where {T <: AbstractFloat} = false
canconvert(::Type{ImagT}, ::Type{Complex{T}}) where {T <: Real} = true
canconvert(::Type{ImagT}, ::Type{Complex}) = true
canconvert(::Type{ImagT}, ::Type{Real}) = false

###
### RootImag
###

canconvert(::Type{RootImagT}, ::Type{T}) where {T <: Real} = false
canconvert(::Type{RootImagT}, ::Type{T}) where {T <: Integer} = false
canconvert(::Type{RootImagT}, ::Type{Rational}) = false
canconvert(::Type{RootImagT}, ::Type{Rational{T}}) where {T <: Integer} = false
canconvert(::Type{RootImagT}, ::Type{T}) where {T <: AbstractFloat} = false
canconvert(::Type{RootImagT}, ::Type{Complex{T}}) where {T <: Rational} = false
canconvert(::Type{RootImagT}, ::Type{Complex{T}}) where {T <: Real} = true
canconvert(::Type{RootImagT}, ::Type{Complex}) = true
canconvert(::Type{RootImagT}, ::Type{Real}) = false

###
### Define the conversions. Rely on `canconvert` above to determine whether to write a method.
###

Base.Bool(::OneT) = true
Base.Bool(::ZeroT) = false

function _make_types()
    intypes = [:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128, :BigInt]
    floatypes = [:Float32, :Float64, :BigFloat]
    rattypes = [:(Rational{$t}) for t in intypes]
    comptypes = [:(Complex{$t}) for t in intypes]
    compfloattypes = [:(Complex{$t}) for t in floatypes]
    comprattypes = [:(Complex{Rational{$t}}) for t in intypes]

    alltypes = [intypes..., floatypes..., rattypes..., comptypes..., comprattypes..., compfloattypes...]
    return alltypes
end

import Base: Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128, BigInt
import Base: Float32, Float64, BigFloat

# const ALLTYPES = _make_types()
# @show ALLTYPES
# TODO: finding bugs because of lack of testing
import Base: sqrt, AbstractFloat, Complex, Rational, Real, Integer, big

# UPDATE: I fixed some of this with the `continue` below.
# Test this.
# TODO: sqrt is applied when getting a numerica value.
# This can produce a float from an int.
# But this happens when requesting an integer.
# Eg, the following function is written, but it is wrong.
# function Complex{Int}(RootImag)
#     sqrt(Complex{Int}(im)
# end

for (ST, litST, func) in ((:OneT, 1, :identity), (:TwoT, 2, :identity), (:InvTwoT, 1//2, :identity), (:RootTwoT, 2, :sqrt),
                          (:InvRootTwoT, 1//2, :sqrt), (:ImagT, im, :identity), (:RootImagT, im, :sqrt))
    STT = @eval $ST
    for T in _make_types()
        TT = @eval $T
        if (TT <: Integer || TT <: Complex{<:Integer}) && func === :sqrt
            # Caller explicitly wants integer or gaussian integer. But sqrt would convert to float.
            # So skip writing the method
            continue
        end
        if canconvert(STT, TT)
            @eval function ($T)(obj::$ST)
                $func($T($litST))
            end
        end
    end
    ## TODO: Test following. It is wrong in some cases
    if canconvert(STT, AbstractFloat)
        @eval function AbstractFloat(obj::$ST)
       #            $func(Float64(obj)) # this appies func twice
                    Float64(obj)
        end
    end
    if canconvert(STT, Rational)
        @eval function Rational(obj::$ST)
                Rational{Int64}(obj)
        end
    end
    if canconvert(STT, Int64)
        @eval function Real(obj::$ST)
            Int64(obj)
        end
        @eval function Integer(obj::$ST)
            Int64(obj)
        end
        @eval function big(obj::$ST)
            BigInt(obj)
        end
    elseif canconvert(STT, Rational{Int64})
        @eval function Real(obj::$ST)
             Rational{Int64}(obj)
        end
        @eval function big(obj::$ST)
            Rational{BigInt}(obj)
        end
    elseif canconvert(STT, Float64)
        @eval function Real(obj::$ST)
             Float64(obj)
        end
        @eval function big(obj::$ST)
            BigFloat(obj)
        end
    end

    notreal = ! canconvert(STT, Real)
    if canconvert(STT, Complex{Int64})
        @eval function Complex(obj::$ST)
                Complex{Int64}(obj)
        end
        if notreal
            @eval function big(obj::$ST)
                Complex{BigInt}(obj)
            end
        end
    elseif canconvert(STT, Complex{Rational{Int64}})
        @eval function Complex(obj::$ST)
                Complex{Rational{Int64}}(obj)
        end
        if notreal
            @eval function big(obj::$ST)
                Complex{Rational{BigInt}}(obj)
            end
        end
    elseif canconvert(STT, Complex{Float64})
        @eval function Complex(obj::$ST)
                Complex{Float64}(obj)
        end
        if notreal
            @eval function big(obj::$ST)
                Complex{BigFloat}(obj)
            end
        end
    end
end

# We could make promote_rule rules.
# We use `Number` to avoid dispatching to these when calling
# with singletons

Base.:*(pow::Pow{T}, x::Pow) where {T} = throw(MethodError(Base.:*, (pow, x)))
Base.:*(pow::Pow{T}, x::NT) where {NT <: Number, T} = NT(T())^pow.n * x
Base.:*(x::Number, pow::Pow) = x * pow

function Base.literal_pow(f::typeof(^), x::SingleNum, ::Val{N}) where {N}
    Pow{typeof(x)}(N)
end

function Base.:^(x::SingleNum, n::Integer)
    Pow{typeof(x)}(n)
end
Base.:^(::OneT, n::Integer) = One

Base.:^(::ImagT, n::Integer) = Pow{ImagT}(mod(n, 4))
function Base.literal_pow(f::typeof(^), x::ImagT, ::Val{N}) where {N}
    Pow{typeof(x)}(mod(N, 4))
end

Base.:^(::RootImagT, n::Integer) = Pow{RootImagT}(mod(n, 8))
function Base.literal_pow(f::typeof(^), x::RootImagT, ::Val{N}) where {N}
    Pow{typeof(x)}(mod(N, 8))
end

Base.:*(x::Pow{TwoT}, ::InvTwoT) = Pow{TwoT}(x.n - 1)
Base.:*(x::Pow{InvTwoT}, ::TwoT) = Pow{InvTwoT}(x.n + 1)
Base.:*(x::Pow{TwoT}, y::Pow{TwoT}) = Pow{TwoT}(x.n + y.n)
Base.:*(x::Pow{TwoT}, y::Pow{InvTwoT}) = Pow{TwoT}(x.n - y.n)
Base.:*(x::Pow{InvTwoT}, y::Pow{TwoT}) = Pow{TwoT}(y.n - x.n)
Base.:*(x::Pow{InvTwoT}, y::Pow{InvTwoT}) = Pow{InvTwoT}(x.n + y.n)

Base.:*(x::Pow{TwoT}, y::Pow{InvRootTwoT}) = Pow{InvRootTwoT}(y.n - 2 * x.n)
Base.:*(x::Pow{InvRootTwoT}, y::Pow{TwoT}) = y * x

Base.:*(x::Pow{InvRootTwoT}, y::Pow{InvRootTwoT}) = Pow{InvRootTwoT}(x.n + y.n)

Base.:*(::InvRootTwoT, y::Pow{InvRootTwoT}) = Pow{InvRootTwoT}(y.n + 1)

Base.inv(x::Pow{T}) where {T} = inv(T()) ^ (x.n)
Base.inv(::OneT) = One

Base.:*(::ZeroT, ::ZeroT) = Zero
Base.:*(::ZeroT, ::SingleNum) = Zero
Base.:*(::SingleNum, ::ZeroT) = Zero
Base.:*(::OneT, ::ZeroT) = Zero
Base.:*(::ZeroT, ::OneT) = Zero
Base.:*(::OneT, s::SingleNum) = s
Base.:*(s::SingleNum, ::OneT) = s
Base.:*(::OneT, ::OneT) = One
Base.:*(::OneT, x::Number) = x
Base.:*(x::Number, ::OneT) = x
Base.:*(::RootImagT, ::RootImagT) = Imag
Base.:*(::RootTwoT, ::RootTwoT) = Two
Base.:*(::InvRootTwoT, ::InvRootTwoT) = InvTwo
# Not using OneT here.
Base.:*(::InvTwoT, ::TwoT) = One
Base.:*(::TwoT, ::InvTwoT) = One
Base.:*(::InvRootTwoT, ::RootTwoT) = One
Base.:*(::RootTwoT, ::InvRootTwoT) = One
Base.:*(::ImagT, ::ImagT) = -1

function Base.:*(::RootImagT, x::T) where {T<:Number}
    T2 = complex(T)
    (T2(zero(T), one(T)) + one(T2)) / sqrt(T2(2)) * x
end

Base.:*(x::SingleNum, y::Number) = typeof(y)(x) * y

Base.:*(::OneT, ::RootImagT) = RootImag
Base.:*(x::SingleNum, ::RootImagT) = RootImag * x
Base.:*(::ZeroT, ::RootImagT) = Zero

Base.:*(::ImagT, x::T) where {T<:Number} = complex(T)(im) * x
# Base.:*(x, ::ImagT) = Imag * x

Base.:+(::OneT, ::OneT) = Two
Base.:-(::TwoT, ::OneT) = One
Base.:-(::T, ::T) where {T <: SingleNum} = Zero
Base.:-(::ZeroT, ::ZeroT) = Zero
Base.:-(x::NumSingleNum, ::ZeroT) = x
Base.:+(::ZeroT, x::NumSingleNum) = x

# Don't know where this is used
Base.promote_rule(::Type{TwoT}, ::Type{T}) where {T<:Integer} = T

# Base.:*(::TwoT, x::T) where {T} = T(2) * x
#Base.:*(x, ::TwoT) = Two * x

Base.:*(::InvRootTwoT, x::T) where {T<:Number} = sqrt(float(T)(1//2)) * x
#Base.:*(x, ::InvRootTwoT) = InvRootTwo * x

end # module SingletonNumbers
