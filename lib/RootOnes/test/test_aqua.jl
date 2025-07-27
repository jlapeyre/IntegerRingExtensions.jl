using RootOnes
using Aqua: Aqua

const PACKAGE_NAME = RootOnes

# TODO:
@testset "aqua test ambiguities $PACKAGE_NAME" begin
    Aqua.test_ambiguities([PACKAGE_NAME, Core, Base])
end

@testset "aqua unbound_args $PACKAGE_NAME" begin
    Aqua.test_unbound_args(PACKAGE_NAME)
end

@testset "aqua undefined exports $PACKAGE_NAME" begin
    Aqua.test_undefined_exports(PACKAGE_NAME)
end

@testset "aqua project extras $PACKAGE_NAME" begin
    Aqua.test_project_extras(PACKAGE_NAME)
end

@testset "aqua stale deps $PACKAGE_NAME" begin
    Aqua.test_stale_deps(PACKAGE_NAME)
end

@testset "aqua deps compat $PACKAGE_NAME" begin
    Aqua.test_deps_compat(PACKAGE_NAME)
end

@testset "aqua piracies $PACKAGE_NAME" begin
    Aqua.test_piracies(PACKAGE_NAME)
end

# Quite slow. Unlikely we violate this in any case
# @testset "aqua persistent tasks $PACKAGE_NAME" begin
#     Aqua.test_persistent_tasks(PACKAGE_NAME)
# end

# TODO:
# @testset "aqua undocumented_names $PACKAGE_NAME" begin
#     Aqua.test_undocumented_names(PACKAGE_NAME)
# end
