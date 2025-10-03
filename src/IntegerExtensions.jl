module IntegerExtensions
import Reexport

using DispatchDoctor: @stable, @unstable

export imaginary, sqrt_imaginary
export DRoot2
export DOmega
export ZRoot2
export ZRootD

include("utils.jl")
Reexport.@reexport import .Utils: subscript, superscript, lobit, small, countmap, random_unitary,
    random_special_unitary, prettylist

include("common.jl")
Reexport.@reexport import .Common: canonical, one_over_root_two, root_two, imaginary, sqrt_imaginary, coeffs, params,
    mul_root_two, mul_one_over_root_two, mul_half, mul_two, conj_root_two, norm_root_two, norm_root_D, isrational, isunit, invstrict,
    isimag

include("rings/singletons.jl")
Reexport.@reexport import .Singletons: SingleNum, RootTwo, Imag, 𝕚, RootImag, InvRootTwo, Two, 𝟚,
    InvTwo, 𝟚⁻¹, 𝟙, One, Zero, 𝟘, canconvert, Pow,
    TwoT, RootTwoT, ImagT, RootImagT, InvRootTwoT, TwoT, InvTwoT, OneT, ZeroT

include("rings/root_one.jl")
Reexport.@reexport import .RootOnes: RootOne, Omega, omega

include("rings/dyadic.jl")
Reexport.@reexport import .Dyadics: Dyadic

include("rings/quadratic_ring.jl")
Reexport.@reexport import .QuadraticRings: QuadraticRing, ZRootD, ZRoot2, DRoot2, QuadraticRing2,
    conj_root_D

include("rings/cyclotomic_rings.jl")
Reexport.@reexport import .CyclotomicRings: CyclotomicRing, DOmega, DOmegaA, ZOmega, least_denominator_exponent,
    rnorm, div_half

# include("secants.jl")
# Reexport.@reexport import .Secants: secant, secant_slope

include("matrices.jl")
Reexport.@reexport import .Matrices2x2:
    AbstractUnitaryNxN,AbstractNormalNxN,
    AbstractMatrix2x2, AbstractMatrixNxN,  AbstractMatrix4x4,
    AbstractUnitary2x2, AbstractSU2,
    MatrixNxN, Matrix4x4, ScaleMatrix2x2, scalematrix,
    Matrix2x2, Vector2,
    GPID, random_diagonal_unitary,
    tracenorm, tracedistance,
    isSU2, random_unitary2x2,
    SU2, SU2B, SU2C, Unitary2x2P, elements, ZRot, zrot,
    get_theta, unitary_u, unitary_t,
    random_ZRot, columns, opnormdistance, alt_random_unitary2x2,
    isantidiag, eigvals_hermitian, antihermitianpart, isantihermitian,
    eigen_hermitian, check_eigen, matrix_func

include("qm_matrices.jl")
Reexport.@reexport import .QMMatrices: random_density_matrix2x2, DensityMatrix2x2,
    ispure, polarization, polarization_norm, depolarize, isdensitymatrix

include("ring_matrices.jl")
Reexport.@reexport import .RingMatrices: compute_phase_factor

include("gates.jl")
Reexport.@reexport import .Gates: Gate1, RZ, X, Y, Z, I2x2, S, T, H, SX, SY, SH,
    commutator

include("compose.jl")
Reexport.@reexport import .Compose: compose, compose_one, compose_scale, get_global_phase,
    correct_global_phase, rotation_error, rotation_error_GPID

include("gate_matrix.jl")
Reexport.@reexport import .GateMatrix: XF64, YF64, ZF64, SF64, TF64,
    CLIFFORD_STR, CLIFFORD_DOMEGA, CLIFFORD_DOMEGA_MAP, find_clifford

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
