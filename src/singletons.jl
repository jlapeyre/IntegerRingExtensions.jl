module Singletons

import Base: show, inv, sqrt
import ..Utils: PRETTY

export RootTwo, InvRootTwo, Imag, 𝕚,  RootImag, Two, 𝟚, InvTwo, 𝟚⁻¹,
    One, Zero, 𝟙, 𝟘

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

"""
    Zero
    𝟘

Represents the additive identity: 𝟘
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

"""
    One
    𝟙

Represents the multiplicative identity: 𝟙
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
sqrt(Two) = RootTwo

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

"""
    RootImag

Represents the principal square root of the imaginary unit: √𝕚
This is also the principal eight root of one.
"""
RootImag

"""
    canconvert(::Type{T}, ::Type{V})::Bool
    canconvert(obj::T, ::Type{V})::Bool

Returns `true` if we can convert `obj::T` with the call `V(obj)`.

The result should be an object of type `<: V`.
"""
function canconvert end

# If called on an object, call again, on type of object
canconvert(obj::SingleNum, ::Type{V}) where {V} = canconvert(typeof(obj), V)

# Ugh. Minor problem. Can't specify that a type must be abstract or not abstract :(
# If ST can be converted to some Rational type, it can be converted to some AbstractFloat type.
# That is: If I can convert with `Rational(obj)`, then I can convert with `AbstractFloat(obj)`.
# function canconvert(::Type{ST}, ::Type{AbstractFloat}) where {ST<:SingleNum}
#     canconvert(ST, Rational)
# end

###
### <: SingleNum
###


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

for ST in (:RootTwoT, :InvRootTwoT)
    @eval canconvert(::Type{$ST}, ::Type{T}) where {T <: Integer} = false
    @eval canconvert(::Type{$ST}, ::Type{T}) where {T <: AbstractFloat} = true
end

###
### InvTwo
###

canconvert(::Type{InvTwoT}, ::Type{T}) where {T <: Integer} = false
canconvert(::Type{InvTwoT}, ::Type{Rational}) = true

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

###
### Maybe most of the stuff below is not necessary.
###

## These are DANGEROUS. Someone may be tempted to use them.
## Used dynamically, they will give terrible performance
## Maybe better to use a more generic function and if-else to distinguish
## elements.

Base.:*(::RootImagT, ::RootImagT) = Imag
Base.:*(::RootTwoT, ::RootTwoT) = Two
Base.:*(::InvRootTwoT, ::InvRootTwoT) = InvTwo

Base.:*(::InvTwoT, ::TwoT) = 1
Base.:*(::TwoT, ::InvTwoT) = 1

Base.:*(::InvRootTwoT, ::RootTwoT) = 1 # May want to use a different type
Base.:*(::RootTwoT, ::InvRootTwoT) = 1 # May want to use a different type

Base.:*(::ImagT, ::ImagT) = -1

function _make_type_pairs()
    _pairs = []
    for s in (:Int8, :Int16, :Int32, :Int64, :Int128, :UInt8, :UInt16, :UInt32, :UInt64, :UInt128)
        push!(_pairs, (s, s))
        rs = :(Rational{$s})
        push!(_pairs, (rs, rs))
        rs = :(Complex{$s})
        push!(_pairs, (rs, rs))
    end
    for s in (:Float64, :Float32, :Float16)
        push!(_pairs, (s, s))
        rs = :(Complex{$s})
        push!(_pairs, (rs, rs))
    end
    for ps in ((:Integer, :Int64), (:Rational, :(Rational{Int})), (:AbstractFloat, :Float64),)
        push!(_pairs, ps)
    end
    _pairs
end

for (Ta, Tb) in _make_type_pairs()
    for (V, val) = ((:TwoT, 2), (:InvTwoT, 1//2), (:ImagT, im))

        tT = @eval $Ta
        if V in (:InvTwoT,) && (tT <: Integer || tT <: Complex{<:Integer})
            continue
        end
        if V in (:ImagT,) && !(tT <: Complex)
            continue
        end
        @eval function Base.convert(::Type{$Ta}, ::$V)
            $Tb($val)
        end
        @eval function (::Type{$Ta})(::$V)
            $Tb($val)
        end
    end
end

function Base.:*(::RootImagT, x::T) where {T}
    T2 = complex(T)
    (T2(im) + one(T2))/sqrt(T2(2)) * x
end

Base.:*(x, ::RootImagT) = RootImag * x

Base.:*(::ImagT, x::T) where {T} = complex(T)(im) * x
Base.:*(x, ::ImagT) = Imag * x

Base.:*(::OneT, x::T) where {T} = one(T) * x
Base.:*(x::T, ::OneT) where {T} = x * one(T)

Base.:*(::InvRootTwoT, x::T) where {T} = sqrt(float(T)(1//2)) * x
Base.:*(x::T, ::InvRootTwoT) where {T} = InvRootTwo * x

end # module Singletons
