using IntegerExtensions
using Test

push!(LOAD_PATH, joinpath(@__DIR__, "..", "lib"))

# Angles2
include(joinpath(@__DIR__, "..", "lib", "Angles2", "test", "runtests.jl"))

# Functionsecants
include(joinpath(@__DIR__, "..", "lib", "FunctionSecants", "test", "runtests.jl"))

include("matrices_tests.jl")
include("rootone_tests.jl")
include("singleton_conversion_tests.jl")
include("composition_tests.jl")
include("quadratic_ring_tests.jl")
include("cyclotomic_tests.jl")
include("dyadic_tests.jl")
include("common_tests.jl")
include("zrot_tests.jl")
include("ringmatrices_tests.jl")
# include("test_aqua.jl")
