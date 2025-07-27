using IntegerExtensions
import LinearAlgebra: norm
using Test

@testset "QuadraticRing" begin
    q = QuadraticRing(1, 1, 2)
    @test q === QuadraticRing{2}(1, 1)
    @test iszero(zero(q))
    @test iszero(zero(QuadraticRing{2, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm(q^5) == norm(q)
end
