@testset "RootOne" begin
    # Test that `angle` gives the same result obtained by converting to
    # float, and then extracting the angle.
    for i in 0:7
        root = RootOne8(i)
        @test isapprox(angle(root), angle(float(root)))
        @test isone(root * inv(root))
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
end
