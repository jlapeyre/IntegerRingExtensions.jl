using IntegerExtensions
using Test

include("singleton_conversion_tests.jl")
include("composition_tests.jl")
include("quadratic_ring_tests.jl")
include("cyclotomic_tests.jl")
include("dyadic_tests.jl")
include("common_tests.jl")
include("matrices_tests.jl")
include("zrot_tests.jl")
include("ringmatrices_tests.jl")
include("rootonetests.jl")

@testset "QuadraticRing{2, Int}" begin
    # Don't follow this example in real code!
    # If D is not literal or `const`, performance degrades by orders of magnitude!
    D = 2
    q = QuadraticRing{D}(1, 1)
    @test q === QuadraticRing{D}(1, 1) # different ways to instantiate.
    @test iszero(zero(q))
    @test iszero(zero(QuadraticRing{D, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm_root_two(q^5) == norm_root_two(q)
    @test conj_root_two(q) * q == norm_root_two(q)

    q2 = QuadraticRing2(1, 2)
    @test norm_root_two(q2) == -7
    @test conj_root_two(q2) * q2 == norm_root_two(q2)
    @test conj_root_two(q2) * q2 == -7.0

    u1 = QuadraticRing{D}(1, 1)
    u2 = QuadraticRing{D}(1, 0)
    @test isunit(u1)  # Units are not unique
    @test isunit(u2)

    @test convert(Int, QuadraticRing{D}(3, 0)) == 3
    @test_throws ArgumentError convert(Int, QuadraticRing{D}(3, 1)) # Can't convert 3 + sqrt(2) to Int
    @test float(QuadraticRing{D}(3, 2)) == 3 + 2 * sqrt(2)
end


@testset "DOmega" begin
    cr = CyclotomicRing(1,2,3,4)
    @test isa(cr, CyclotomicRing{4, Int64})
    cromega = DOmega((1,2,3,4))
    @test cromega === DOmega(1,2,3,4)
    @test typeof(cromega) === DOmega{Int}
    @test cromega === DOmega{Int}(cr)
    z = DOmega(1,2,3,Dyadic(1,1))
    @test typeof(z) === DOmega{Int}
    @test typeof(DOmega(1,2,3,4)) === DOmega{Int}

    @test cromega == cromega
    @test DOmega(1,0,0,0) == one(cromega)
    @test isone(one(cr))
    @test isone(one(cromega))

    azf = float(z * conj(z))
    @test isreal(azf)
    @test real(azf) > 0
end

@testset "DRoot2" begin
    x = DRoot2(1, Dyadic(3, 2))
    @test typeof(x) === DRoot2{Int64, Int64}
    @test typeof(big(x)) === QuadraticRing{2, BigFloat}
    @test typeof(big(big(x))) === BigFloat
end

@testset "Matrix{QuadraticRing}" begin
    # z = 1/sqrt(2)
    z = QuadraticRing{2}(Dyadic(0,0), Dyadic(1,1))
    @test float(z) == sqrt(2) / 2
    hadamard = [z z; z -z] # Hadamard matrix
    @test isapprox(float(hadamard), 1/sqrt(2) * [1 1; 1 -1])
end
