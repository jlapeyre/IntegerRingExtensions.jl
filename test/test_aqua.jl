using IntegerExtensions
using Aqua: Aqua

# TODO:
# @testset "aqua test ambiguities" begin
#     Aqua.test_ambiguities([IntegerExtensions, Core, Base])
# end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(IntegerExtensions)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(IntegerExtensions)
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(IntegerExtensions)
end

@testset "aqua stale deps" begin
    Aqua.test_stale_deps(IntegerExtensions)
end

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(IntegerExtensions)
end

@testset "aqua piracies" begin
    Aqua.test_piracies(IntegerExtensions)
end

@testset "aqua persistent tasks" begin
    Aqua.test_persistent_tasks(IntegerExtensions)
end

# TODO:
# @testset "aqua undocumented_names" begin
#     Aqua.test_undocumented_names(IntegerExtensions)
# end
