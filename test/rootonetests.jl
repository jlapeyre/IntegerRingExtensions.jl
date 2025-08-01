@testset "RootOne" begin
    # Test that `angle` gives the same result obtained by converting to
    # float, and then extracting the angle.
    for i in 0:7
        root = RootOne8(i)
        @test isapprox(angle(root), angle(float(root)))
        @test isone(root * inv(root))
        @test isone(inv(root) * root)
        @test isone(conj(root) * root)
        @test conj(root) == adjoint(root)
    end
    @test RootOne8 === RootOne{8}
    @test RootOne8() === RootOne8(1)
    r = RootOne8(1)
    @test isa(angle(r), Float64)
    @test isa(angle(BigFloat, r), BigFloat)

    # Test that BigFloat angle has full precision
    r5 = RootOne8(5)
    z = cis(angle(BigFloat, r5))
    @test real(z^2) < 1 / big(10)^75

    @test RootOne8(4) / RootOne8(2) === RootOne8(2)

    @test RootOne{8}(2) == RootOne{4}(1)
    @test RootOne{5}(0) == RootOne{17}(0)

    @test imaginary(RootOne{4}) == RootOne{4}(1)
    @test imaginary(RootOne{8}) == RootOne{8}(2)
    @test imaginary(RootOne{12}) == RootOne{12}(3)
    @test_throws ArgumentError imaginary(RootOne{6})

    @test sqrt_imaginary(RootOne{8}) == RootOne{8}(1)
    @test sqrt_imaginary(RootOne{16}) == RootOne{16}(2)
    @test_throws ArgumentError sqrt_imaginary(RootOne{12})
end
