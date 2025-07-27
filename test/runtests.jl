using IntegerRingExtensions
using Test

include("test_aqua.jl")

push!(LOAD_PATH, joinpath(@__DIR__, "..", "lib"))

# Angles2
include(joinpath(@__DIR__, "..", "lib", "Angles2", "test", "runtests.jl"))

# No tests here
# FunctionSecants
#include(joinpath(@__DIR__, "..", "lib", "FunctionSecants", "test", "runtests.jl"))

# SingletonNumbers
include(joinpath(@__DIR__, "..", "lib", "SingletonNumbers", "test", "runtests.jl"))

# RootOnes
include(joinpath(@__DIR__, "..", "lib", "RootOnes", "test", "runtests.jl"))

# Dyadics
include(joinpath(@__DIR__, "..", "lib", "Dyadics", "test", "runtests.jl"))

# QuadraticRings
include(joinpath(@__DIR__, "..", "lib", "QuadraticRings", "test", "runtests.jl"))

# CyclotomicRings
include(joinpath(@__DIR__, "..", "lib", "CyclotomicRings", "test", "runtests.jl"))

include("matrices_tests.jl")
include("composition_tests.jl")
include("common_tests.jl")
include("zrot_tests.jl")
include("ringmatrices_tests.jl")
