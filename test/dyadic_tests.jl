@testset "Dyadics" begin
    @test Dyadic(4, 2) == Dyadic(1, 0)
    @test Int(Dyadic(4, 1)) == 2
    @test Int(Dyadic(4, 0)) == 4
    @test Int(Dyadic(4, 2)) == 1
    @test Rational(Dyadic(3, 3)) === 3//8
    @test Dyadic(3//8) === Dyadic(3, 3)
    @test_throws ArgumentError Int(Dyadic(3, 2))

    z = Dyadic(3, 1)
    @test float(z) == 3/2

    zh = mul_half(z)
    @test zh === Dyadic(3, 2)
    @test float(zh) == 3/4

    z = Dyadic(4, 1)
    zh = mul_half(z)
    @test float(zh) == float(z) / 2
    cs = params(zh)
    @test cs == (1, 0)

    x = Droot2(Dyadic(3,2), 5)
    @test RootTwo * x == Droot2(10, Dyadic(3,2))
    @test x * RootTwo == Droot2(10, Dyadic(3,2))

    @test x === InvRootTwo * (RootTwo * x)
    @test x === RootTwo * (InvRootTwo * x)

    @test isinteger(Dyadic(3, 0))
    @test !isinteger(Dyadic(3, 1))
end
