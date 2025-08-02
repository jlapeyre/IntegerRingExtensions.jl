module Singletons

import Base: show, inv, sqrt
import ..Utils: PRETTY

export RootTwo, InvRootTwo, Imag, 𝕚,  RootImag, Two, 𝟚, InvTwo, 𝟚⁻¹,
    𝟙, One, Zero, 𝟘

###
### Zero
###

struct ZeroT end
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

struct OneT end
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

struct TwoT end
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

struct InvTwoT end
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

struct RootTwoT end
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

struct InvRootTwoT end
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

struct ImagT end
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

struct RootImagT end
const RootImag = RootImagT()
show(io::IO, ::PRETTY, ::RootImagT) = print(io, "√𝕚")
sqrt(::ImagT) = RootImag

"""
    RootImag

Represents the principal square root of the imaginary unit: √𝕚
This is also the principal eight root of one.
"""
RootImag

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

Base.:*(::ImagT, x::T) where {T} = T(im) * x
Base.:*(x, ::ImagT) = Imag * x

Base.:*(::OneT, x::T) where {T} = one(T) * x
Base.:*(x::T, ::OneT) where {T} = x * one(T)

end # module Singleons
