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

    z = DyadicFraction(4, 1)
    zh = mul_half(z)
    @test float(zh) == float(z) / 2
    cs = params(zh)
    @test cs == (1, 0)

    x = Droot2(DyadicFraction(3,2), 5)
    @test RootTwo * x == Droot2(10, DyadicFraction(3,2))
    @test x * RootTwo == Droot2(10, DyadicFraction(3,2))

    @test x === InvRootTwo * (RootTwo * x)
    @test x === RootTwo * (InvRootTwo * x)
end
