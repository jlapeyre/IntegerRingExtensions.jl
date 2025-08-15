module IntegerExtensions
import Reexport

export imaginary, sqrt_imaginary, Droot2
export Domega
export Droot2
export Zroot2
export ZrootD

include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript, lobit, random_angle, small, countmap, random_unitary,
    random_special_unitary, random_phase

include("common.jl")
Reexport.@reexport import .Common: canonical, one_over_root_two, root_two, imaginary, sqrt_imaginary, coeffs, params,
    mul_root_two, mul_one_over_root_two, mul_half, conj_root_two, norm_root_two, isrational, isunit, invstrict

include("singletons.jl")
Reexport.@reexport import .Singletons: SingleNum, RootTwo, Imag, 𝕚, RootImag, InvRootTwo, Two, 𝟚,
    InvTwo, 𝟚⁻¹, 𝟙, One, Zero, 𝟘, canconvert, Pow,
    TwoT, RootTwoT, ImagT, RootImagT, InvRootTwoT, TwoT, InvTwoT, OneT, ZeroT

include("angles.jl")
Reexport.@reexport import .Angles: Dar, dartorad, radtodar, scalepi, unscalepi, minus_one_to_one

include("rings/root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, Omega

include("rings/dyadic.jl")
Reexport.@reexport import .Dyadics: Dyadic

include("rings/quadratic_ring.jl")
Reexport.@reexport import .QuadraticRings: QuadraticRing, ZrootD, Zroot2, Droot2, QuadraticRing2,
    rootDconj

include("rings/cyclotomic_rings.jl")
Reexport.@reexport import .CyclotomicRings: CyclotomicRing, Domega, Zomega, smallest_denominator_exponent

include("matrices.jl")
Reexport.@reexport import .Matrices2x2: Matrix2x2, Vector2, GPID, random_diagonal_unitary,
    tracenorm, tracedistance, unitary_compose, unitary_decompose,
    UnitaryParam1, UnitaryParam2,
     isSU2, random_unitary2x2, random_SU2,
    SU2OLD, Unitary2x2, elements, ZRot, zrot,
    get_theta, SU2, unitary_u, unitary_t,
    AbstractUnitary2x2, AbstractSU2, AbstractMatrix2x2,
    random_ZRot
#    zrotpi, zrothalfpi, get_thetapi, get_thetahalfpi
#    SU2Param1, SU2Param2, SU2Param3, SU2ParamScaled,

include("ring_matrices.jl")
Reexport.@reexport import .RingMatrices: compute_phase_factor

include("compose.jl")
Reexport.@reexport import .Compose: compose, compose_one, RZ, Gate1, get_global_phase,
    correct_global_phase, rotation_error, rotation_error_GPID, Uapprox, random_Uapprox, random_RZ,
    isUapprox

include("benchmark.jl")
Reexport.@reexport import .Benchmarking: benchmark_compose

include("gridsynth.jl")
Reexport.@reexport import .GridSynth: GridSynthOpts, GridSynthResults, makecommand, run_gridsynth, gridsynth,
    stringtonum, gridsynth_matrix, biggenex, biggennum

include("gridsynth_extra.jl")
Reexport.@reexport import .GridSynthExtra: gridsynth_matrix_to_cyclic

# For compiling workflows for statically-compiled-like latency
# using PrecompileTools: @setup_workload, @compile_workload
# include("precompile.jl")

end # module IntegerExtensions
