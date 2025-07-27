using CyclotomicRings
using Aqua: Aqua

const PACKAGE_NAME = CyclotomicRings

# TODO:
# @testset "aqua test ambiguities" begin
#     Aqua.test_ambiguities([PACKAGE_NAME, Core, Base])
# end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(PACKAGE_NAME)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(PACKAGE_NAME)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(PACKAGE_NAME)
end

@testset "aqua stale deps" begin
    Aqua.test_stale_deps(PACKAGE_NAME)
end

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(PACKAGE_NAME)
end

# We intentionally pirate methods that in principle belong
# to RootOnes. Maybe there is a way to suppress these errors.
# @testset "aqua piracies" begin
#     Aqua.test_piracies(PACKAGE_NAME)
# end

# Quite slow. Unlikely we violate this in any case
# @testset "aqua persistent tasks" begin
#     Aqua.test_persistent_tasks(PACKAGE_NAME)
# end

# TODO:
# @testset "aqua undocumented_names" begin
#     Aqua.test_undocumented_names(PACKAGE_NAME)
# end
