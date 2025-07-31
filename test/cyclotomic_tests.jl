import IntegerExtensions.CyclotomicRings: mul_sqrt2

@testset "cyclotomic" begin
    d = Domega(1,2,3,4)
    z = Zomega(1,2,3,4)
    for i in 0:7
        r = RootOne8(i)
        romega = Domega{Int}(r)
        @test r * d == romega * d
        @test r * z == romega * z
    end
    @test isapprox(float(mul_sqrt2(d)) / float(d), sqrt(2))
    @test isapprox(float(mul_sqrt2(z)) / float(z), sqrt(2))
end
