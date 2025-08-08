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
Reexport.@reexport import .Singletons: SingleNum, RootTwo, Imag, 𝕚, RootImag, InvRootTwo, Two, 𝟚,
    InvTwo, 𝟚⁻¹, 𝟙, One, Zero, 𝟘, canconvert, Pow,
    TwoT, RootTwoT, ImagT, RootImagT, InvRootTwoT, TwoT, InvTwoT, OneT, ZeroT

include("rings/root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, RootOne8

include("rings/dyadic.jl")
Reexport.@reexport import .Dyadics: Dyadic

include("rings/quadratic_ring.jl")
Reexport.@reexport import .QuadraticRings: QuadraticRing, ZrootD, Zroot2, Droot2, QuadraticRing2, isunit,
    rootDconj

include("rings/cyclotomic_rings.jl")
Reexport.@reexport import .CyclotomicRings: CyclotomicRing, Domega, Zomega

include("matrices.jl")
Reexport.@reexport import .Matrices2x2: Matrix2x2, Vector2

include("ring_matrices.jl")
Reexport.@reexport import .RingMatrices: compute_phase_factor

include("gates.jl")
Reexport.@reexport import .Gates: compose, compose_one, RZ, count_gates, Gate1, get_theta, get_global_phase,
    correct_global_phase, rotation_error

include("benchmark.jl")
Reexport.@reexport import .Benchmarking: benchmark_compose

include("gridsynth.jl")
Reexport.@reexport import .GridSynth: GridSynthOpts, GridSynthResults, makecommand, run_gridsynth, gridsynth,
    stringtonum, gridsynth_matrix

include("gridsynth_extra.jl")
Reexport.@reexport import .GridSynthExtra: gridsynth_matrix_to_cyclic

# For compiling workflows for statically-compiled-like latency
# using PrecompileTools: @setup_workload, @compile_workload
# include("precompile.jl")

end # module IntegerExtensions
