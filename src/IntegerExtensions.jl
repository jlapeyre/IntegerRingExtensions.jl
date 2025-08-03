module IntegerExtensions
import Reexport

export imaginary, sqrt_imaginary, Droot2
export Domega
export Droot2
export Zroot2
export ZrootD

include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript, lobit

include("common.jl")
Reexport.@reexport import .Common: canonical, one_over_root_two, root_two, imaginary, sqrt_imaginary, coeffs, params,
    mul_root_two, mul_one_over_root_two, mul_half, conj_root_two, norm_root_two, isrational

include("singletons.jl")
Reexport.@reexport import .Singletons: RootTwo, Imag, 𝕚, RootImag, InvRootTwo, Two, 𝟚,
    InvTwo, 𝟚⁻¹, 𝟙, One, Zero, 𝟘, canconvert

include("root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, RootOne8

include("dyadic.jl")
Reexport.@reexport import .DyadicFractions: DyadicFraction

include("quadratic_ring.jl")
Reexport.@reexport import .QuadraticRings: QuadraticRing, ZrootD, Zroot2, Droot2, QuadraticRing2, isunit,
    rootDconj

include("cyclotomic_rings.jl")
Reexport.@reexport import .CyclotomicRings: CyclotomicRing, Domega, Zomega

include("matrices.jl")
Reexport.@reexport import .Matrices: Matrix2x2, get_theta

include("gates.jl")
Reexport.@reexport import .Gates: compose, compose_one, RZ, count_gates, Gate1

include("benchmark.jl")
Reexport.@reexport import .Benchmarking: benchmark_compose


end # module IntegerExtensions
