using IntegerExtensions
using Test

include("ringmatrices_tests.jl")
include("composition_tests.jl")
include("singleton_conversion_tests.jl")
include("rootonetests.jl")
include("matricestests.jl")
include("dyadic_tests.jl")
include("d_z_root2_tests.jl")
include("cyclotomic_tests.jl")

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


@testset "Domega" begin
    cr = CyclotomicRing(1,2,3,4)
    @test typeof(cr) === CyclotomicRing{4, Int64}
    cromega = Domega{Int}((1,2,3,4))
    @test cromega === Domega{Int}(1,2,3,4)
    @test typeof(cromega) === Domega{Int}
    @test cromega === Domega{Int}(cr)
    z = Domega{Int}(1,2,3,DyadicFraction(1,1))
    @test typeof(z) === Domega{Int}
    @test typeof(Domega(1,2,3,4)) === Domega{Int}

    @test cromega == cromega
#    @test Domega{Int}(0,0,0,1) == one(cromega)
    @test Domega{Int}(1,0,0,0) == one(cromega)
    @test isone(one(cr))
    @test isone(one(cromega))

    azf = float(z * conj(z))
    @test isreal(azf)
    @test real(azf) > 0
end

@testset "Droot2" begin
    x = Droot2(1, DyadicFraction(3, 2))
    @test typeof(x) === Droot2{Int64, Int64}
    @test typeof(big(x)) === QuadraticRing{2, BigFloat}
    @test typeof(big(big(x))) === BigFloat
end

@testset "Matrix{QuadraticRing}" begin
    # z = 1/sqrt(2)
    z = QuadraticRing{2}(DyadicFraction(0,0), DyadicFraction(1,1))
    @test float(z) == sqrt(2) / 2
    hadamard = [z z; z -z] # Hadamard matrix
    @test isapprox(float(hadamard), 1/sqrt(2) * [1 1; 1 -1])
end

