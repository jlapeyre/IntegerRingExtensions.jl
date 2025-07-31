using IntegerExtensions
import LinearAlgebra: norm
using Test

include("cyclotomic_tests.jl")
include("matricestests.jl")

@testset "QuadraticRing{2, Int}" begin
    # Don't follow this example in real code!
    # If D is not literal or `const`, performance degrades by orders of magnitude!
    D = 2
    q = QuadraticRing(1, 1, D)
    @test q === QuadraticRing{D}(1, 1) # different ways to instantiate.
    @test iszero(zero(q))
    @test iszero(zero(QuadraticRing{D, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm(q^5) == norm(q)
    @test root2conj(q) * q == norm(q)

    q2 = QuadraticRing2(1, 2)
    @test norm(q2) == -7
    @test root2conj(q2) * q2 == norm(q2)
    @test root2conj(q2) * q2 == -7.0

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
    @test Domega{Int}(0,0,0,1) == one(cromega)
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

@testset "DyadicFraction" begin
    @test DyadicFraction(4, 2) == DyadicFraction(1, 0)
    @test Int(DyadicFraction(4, 1)) == 2
    @test Int(DyadicFraction(4, 0)) == 4
    @test Int(DyadicFraction(4, 2)) == 1
    @test Rational(DyadicFraction(3, 3)) === 3//8
    @test_throws ArgumentError Int(DyadicFraction(3, 2))
end

@testset "composition" begin

    # Approximation of RZ(pi/16) up to global phase, for epsilon = 1e-20

    str = "HTSHTSHTHTHTHTSHTHTSHTHTSHTSHTSHTSHTSHTSHTHTHTSHTSHTSHTHTSHTSHTHTHTHTHTSHTSHTSHTHTHTSHTHTSHTSHTSHTSHTHTHTHTHTSHTHTSHTHTSHTHTSHTHTHTHTSHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTSHTHTHTHTHTHTSHTHTHTHTHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTHTSHTSHTHTSHTHTSHTSHTHTSHTSHTSHTSHTHTHTHTHTHTSHTHTSHTSHTHTHTSHTHTHTSHTHTSHTSHTSHTHTSHTSHTHTSHTSHTHTHTSHTHTHTSHTSHTSHTHTSHTHTHTSHTSHTHTSHTHTHTSHTHTSHTSHTSHTHTHTSHTHTHTHTHTHTHTHTHTHTSHTSHTHTSHTHTSHTHTSHTSHTSHTSHTSHTHTSHTHTHTSHTSHTHTXS"

    # Use ordinary Int64 as the base type.
    m = compose(str)
    m_float = big(m) # Convert from Domega{Int} to BigFloat
    theta = get_theta(m_float) # Extract theta. There will be a global phase
    theta_diff = abs(theta - big(pi) / 16) # Expected theta is pi/16
    @test theta_diff < 1e-19  # Test accuracy
    @test theta_diff > 1e-20  # Sanity check. We are testing something.
end

include("rootonetests.jl")
