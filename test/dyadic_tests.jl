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

    x = DRoot2(Dyadic(3,2), 5)
    @test RootTwo * x == DRoot2(10, Dyadic(3,2))
    @test x * RootTwo == DRoot2(10, Dyadic(3,2))

    @test x === InvRootTwo * (RootTwo * x)
    @test x === RootTwo * (InvRootTwo * x)

    @test isinteger(Dyadic(3, 0))
    @test !isinteger(Dyadic(3, 1))

    @test !isunit(Dyadic(3//1))
    @test !isunit(Dyadic(5, 3))

    for d in (Dyadic(1, 0), Dyadic(2, 0), Dyadic(2, 1), Dyadic(2, 2),
              Dyadic(32, 0), Dyadic(32, 1), Dyadic(32, 5))
        @test isunit(d)
        @test isone(invstrict(d) * d)
    end

    for d in (Dyadic(3//1), Dyadic(5, 3))
        @test !isunit(d)
        @test_throws ArgumentError invstrict(d)
    end

    @test ispow2(Dyadic(4, 0))
    @test ispow2(Dyadic(4, 1))
    @test ispow2(Dyadic(4, 2))
    @test !ispow2(Dyadic(4, 3))
    @test !ispow2(Dyadic(5, 0))
    @test !ispow2(Dyadic(5, 1))

    d = Dyadic{Int64, Int64}(2, 0) + Dyadic{Int64, Int64}(0, 0)*im
    @test isunit(d)
    @test isone(invstrict(d) * d)
end
