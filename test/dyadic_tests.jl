@testset "DyadicFractions" begin
    @test DyadicFraction(4, 2) == DyadicFraction(1, 0)
    @test Int(DyadicFraction(4, 1)) == 2
    @test Int(DyadicFraction(4, 0)) == 4
    @test Int(DyadicFraction(4, 2)) == 1
    @test Rational(DyadicFraction(3, 3)) === 3//8
    @test_throws ArgumentError Int(DyadicFraction(3, 2))

    z = DyadicFraction(3, 1)
    @test float(z) == 3/2

    zh = mul_half(z)
    @test zh === DyadicFraction(3, 2)
    @test float(zh) == 3/4
end
