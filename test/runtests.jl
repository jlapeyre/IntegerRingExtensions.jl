using IntegerExtensions
import LinearAlgebra: norm
using Test

@testset "QuadraticInteger" begin
    q = QuadraticInteger(1, 1, 2)
    @test q === QuadraticInteger{2}(1, 1)
    @test iszero(zero(q))
    @test iszero(zero(QuadraticInteger{2, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm(q^5) == norm(q)
end
