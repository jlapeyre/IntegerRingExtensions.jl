@testset "RingMatrices" begin
    m = DOmega(1, 2, 3, 4)
    for n in 0:7
        @test inv(RootOne{8}(n)) * (RootOne{8}(n) * m ) == m
    end
end

@testset "Matrix{QuadraticRing}" begin
    z = QuadraticRing{2}(Dyadic(0,0), Dyadic(1,1))
    @test float(z) == sqrt(2) / 2
    hadamard = [z z; z -z] # Hadamard matrix
    @test isapprox(float(hadamard), 1/sqrt(2) * [1 1; 1 -1])
end
