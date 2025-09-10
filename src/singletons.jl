module Singletons

import Base: show, inv, sqrt, isone, iszero, isinteger, iseven, isreal
import ..Utils: PRETTY, superscript
import ..Common: isrational


export RootTwo, InvRootTwo, Imag, 𝕚,  RootImag, Two, 𝟚, InvTwo, 𝟚⁻¹,
    One, Zero, 𝟙, 𝟘

export Pow

###
### Zero
###

# Should this be a Number? What is a Number
"""
    abstract type SingleNum

An abstract type for reprenting single numbers as singleton types.
"""
abstract type SingleNum end

struct ZeroT <: SingleNum
end
const Zero = ZeroT()
const 𝟘 = Zero
show(io::IO, ::PRETTY, ::ZeroT) = print(io, "𝟘")

isrational(::ZeroT) = true
isinteger(::ZeroT) = true
iszero(::ZeroT) = true
isone(::ZeroT) = false
iseven(::ZeroT) = true
isreal(::ZeroT) = true

"""
    Zero
    𝟘

Represents the additive identity, 𝟘, in the ring of integers.

`Zero` also represents the additive identity in any ring or field that
extends the integers.
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
    √𝟚

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
sqrt(::InvTwoT) = InvRoot2
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

"""
    struct Pow{T}
        n::Int
    end


Represents and instance of type `T` raised to the power `n`.

This is intended for singleton types `T`.
"""
struct Pow{T} <: Number
    n::Int
end

function show(io::IO, ::PRETTY, p::Pow{T}) where {T}
    show(io, PRETTY(), T())
    print(io, superscript(p.n))
end

function Base.AbstractFloat(p::Pow{T}) where {T}
    AbstractFloat(T())^p.n
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

Base.convert(::Type{T}, obj::SingleNum) where {T} = T(obj)

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
function canconvert(::Type{ST}, ::Type{T}) where {ST<:SingleNum, T<:Bool}
    false
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
### Two
###

canconvert(::Type{TwoT}, ::Type{T}) where {T <: Integer} = true
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
import Base: sqrt, AbstractFloat, Complex, Rational, Real, Integer, big

for (ST, litST, func) in ((:TwoT, 2, :identity), (:InvTwoT, 1//2, :identity), (:RootTwoT, 2, :sqrt),
                          (:InvRootTwoT, 1//2, :sqrt), (:ImagT, im, :identity), (:RootImagT, im, :sqrt))
    STT = @eval $ST
    for T in _make_types()
        TT = @eval $T
        if canconvert(STT, TT)
            @eval function ($T)(obj::$ST)
                $func($T($litST))
            end
        end
    end
    if canconvert(STT, AbstractFloat)
        @eval function AbstractFloat(obj::$ST)
                $func(Float64(obj))
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

Base.:*(pow::Pow{T}, x::NT) where {NT <: Number, T} = NT(T())^pow.n * x
Base.:*(x::Number, pow::Pow) = x * pow
Base.:^(x::SingleNum, n::Integer) = Pow{typeof(x)}(n)
Base.:^(::OneT, n::Integer) = One
Base.inv(::OneT) = One
Base.:*(::ZeroT, ::SingleNum) = Zero
Base.:*(::SingleNum, ::ZeroT) = Zero
Base.:*(::OneT, ::ZeroT) = Zero
Base.:*(::ZeroT, ::OneT) = Zero
Base.:*(::OneT, s::SingleNum) = s
Base.:*(s::SingleNum, ::OneT) = s

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
    (T2(im) + one(T2))/sqrt(T2(2)) * x
end

Base.:*(x, ::RootImagT) = RootImag * x

Base.:*(::ImagT, x::T) where {T<:Number} = complex(T)(im) * x
# Base.:*(x, ::ImagT) = Imag * x

# Don't know where this is used
Base.promote_rule(::Type{TwoT}, ::Type{T}) where {T<:Integer} = T

# Base.:*(::TwoT, x::T) where {T} = T(2) * x
#Base.:*(x, ::TwoT) = Two * x

Base.:*(::InvRootTwoT, x::T) where {T<:Number} = sqrt(float(T)(1//2)) * x
#Base.:*(x, ::InvRootTwoT) = InvRootTwo * x

end # module Singletons
