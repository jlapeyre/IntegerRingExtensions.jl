using IntegerExtensions
import LinearAlgebra: norm
using Test

@testset "QuadraticRing{2, Int}" begin
    D = 2
    q = QuadraticRing(1, 1, D)
    @test q === QuadraticRing{D}(1, 1)
    @test iszero(zero(q))
    @test iszero(zero(QuadraticRing{D, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm(q^5) == norm(q)

    u1 = QuadraticRing{D}(1, 1)
    u2 = QuadraticRing{D}(1, 1)
    @test isunit(u1)
    @test isunit(u2)

    @test convert(Int, QuadraticRing{D}(3, 0)) == 3
    @test_throws ArgumentError convert(Int, QuadraticRing{D}(3, 1))
end

@testset "RootOne" begin
    for i in 0:7
        @test isapprox(angle(RootOne8(i)), angle(float(RootOne8(i))))
    end
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
end
