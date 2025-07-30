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

include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript

include("common.jl")
Reexport.@reexport import .Common: canonical

include("singletons.jl")
Reexport.@reexport import .Singletons: Root2, RootI, InvRoot2

include("root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, RootOne8

include("dyadic.jl")
Reexport.@reexport import .DyadicFractions: DyadicFraction

include("rings.jl")
Reexport.@reexport import .Rings: QuadraticRing, isunit, CyclotomicRing,
    QuadraticRing2, Domega, Droot2, root2conj, rootDconj, ZrootD

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, get_theta

include("gates.jl")
Reexport.@reexport import .Gates: Igate, Zgate, Sgate, Tgate, Hgate, Xgate, Wgate, compose, compose_one, gate_map,
    GATE_MAP_BIG_INT, GATE_MAP_INT, GATE_MAP_ZZ, GATE_MAP_BIG_FLOAT, GATE_MAP_INT128,
    RZ

end # module IntegerExtensions
