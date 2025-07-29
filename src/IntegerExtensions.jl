module IntegerExtensions
import Reexport

export Domega, imaginary, sqrt_imaginary, one_over_root_two

"""
    imaginary(::Type{T}) where {T <: Real}
    imaginary(::Type{Complex{T}}) where T

The imaginary unit as type `T`.

If the imaginary unit cannot be represented as `T`, for example
when `T` is `Real`, then `Complex{T}` is returned.
"""
function imaginary(::Type{T}) where {T <: Real}
    complex(zero(T), one(T))
end

function imaginary(::Type{Complex{T}}) where T
    imaginary(T)
end

"""
    sqrt_imaginary(::Type{T})

The principal square root of the imaginary unit as type `T`.

If the imaginary unit cannot be represented as `T`, for example
when `T` is `Real`, then `Complex{T}` is returned.
"""
function sqrt_imaginary(::Type{T}) where {T <: Real}
    cispi(one(T)/4)
end

function sqrt_imaginary(::Type{Complex{T}}) where {T <: Real}
    cispi(one(T)/4)
end

function sqrt_imaginary(::Type{BigFloat})
    cispi(one(BigFloat)/4)
end

function one_over_root_two(::Type{BigFloat})
    one(BigFloat) / sqrt(big(2))
end

canonical(x) = x

include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript

include("rings.jl")
Reexport.@reexport import .Rings: QuadraticRing, isunit, RootOne, DyadicFraction, CyclotomicRing, RootOneA,
    RootOne8, QuadraticRing2, canonical, Domega

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, get_theta

include("gates.jl")
Reexport.@reexport import .Gates: Igate, Zgate, Sgate, Tgate, Hgate, Xgate, compose, gate_map,
    GATE_MAP_BIG_INT, GATE_MAP_INT, GATE_MAP_ZZ, GATE_MAP_BIG_FLOAT

end # module IntegerExtensions
