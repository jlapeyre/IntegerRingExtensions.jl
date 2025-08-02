import IntegerExtensions.CyclotomicRings: mul_root_two

@testset "cyclotomic" begin
    d = Domega(1,2,3,4)
    z = Zomega(1,2,3,4)
    for i in 0:7
        r = RootOne8(i)
        romega = Domega{Int}(r)
        @test r * d == romega * d
        @test r * z == romega * z
    end
    @test isapprox(float(mul_root_two(d)) / float(d), sqrt(2))
    @test isapprox(float(mul_root_two(z)) / float(z), sqrt(2))

end

@testset "cyclotomic other" begin
    z5 = CyclotomicRing{5}(0,1,0,0,0)
    @test float(z5) == float(RootOne{10}(1))

    z5 = CyclotomicRing{5}(0,0,1,0,0)
    @test float(z5) == float(RootOne{10}(2))

    z5 = CyclotomicRing{5}(0,0,0,1,0)
    @test float(z5) == float(RootOne{10}(3))

    z5 = CyclotomicRing{5}(0,0,0,0,1)
    @test float(z5) == float(RootOne{10}(4))
end

@testset "cyclotomic storage" begin
    z = Domega(1,2,3,DyadicFraction(5,2))
    @test isbits(z)
    @test isa(z, Domega{Int})
    zb = CyclotomicRing{4}(1,2,3,DyadicFraction(big(5),2))
    @test !isbits(zb)
    @test isa(zb, Domega{BigInt})

    # Ugh. Not really coerced to a Domega.
    # Second Int type is also Int128 upon promotion
    z2 = Domega(Int128(1),2,3,DyadicFraction(5,2))
    @test isbits(z2)
    @test isa(z2, CyclotomicRing{4, DyadicFraction{Int128,Int128}})
end

@testset "CyclotomicRing construction" begin
    @test Zomega(3) isa Zomega{Int}
    @test Domega(10) isa Domega{Int}
    @test Domega{Int32}(10) isa Domega{Int32}
    # These test both construction and conversion
    @test Int(Zomega(3)) === 3
    @test Int(Domega(3)) === 3

    @test_throws MethodError Domega(1,2,3)
    @test_throws MethodError Zomega(1,2,3)
end
