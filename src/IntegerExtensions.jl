module IntegerExtensions
import Reexport

export imaginary, sqrt_imaginary, Droot2
export Domega
export Droot2
export Zroot2
export ZrootD


include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript

include("common.jl")
Reexport.@reexport import .Common: canonical, one_over_root_two, root_two, imaginary, sqrt_imaginary, coeffs, params,
    mul_root_two, mul_one_over_root_two, mul_half

include("singletons.jl")
Reexport.@reexport import .Singletons: Root2, RootI, InvRoot2

include("root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, RootOne8

include("dyadic.jl")
Reexport.@reexport import .DyadicFractions: DyadicFraction

include("quadratic_ring.jl")
Reexport.@reexport import .QuadraticRings: QuadraticRing, ZrootD, Zroot2, Droot2, QuadraticRing2, isunit,
    rootDconj, conj_root_two

include("cyclotomic_rings.jl")
Reexport.@reexport import .CyclotomicRings: CyclotomicRing, Domega, Zomega

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, get_theta

include("gates.jl")
Reexport.@reexport import .Gates: Igate, Zgate, Sgate, Tgate, Hgate, Xgate, Ygate, Wgate, compose, compose_one, gate_map,
    GATE_MAP_BIG_INT, GATE_MAP_INT, GATE_MAP_BIG_FLOAT, GATE_MAP_INT128, RZ,
    count_gates

end # module IntegerExtensions
