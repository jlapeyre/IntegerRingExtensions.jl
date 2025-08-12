using IsApprox: isunitary

@testset "ZRot" begin
    for theta in (1/2, 1/4, 0.5678, Dar(1/4), Dar(1//4))
        rz = zrot(theta)
        @test isapprox(theta, get_theta(rz))
        ang = - 2 * angle(rz[1])
        @test isapprox(ang, get_theta(rz))
        @test isunitary(rz)
        @test isSU2(rz)
    end

    rz1 = zrot(Dar(1//8))
    rz2 = zrot(Dar(3//8))
    rz3 = rz1 * rz2
    @test get_theta(rz3) === Dar(1//2)

    irz1 = inv(rz1)
    @test typeof(irz1) == typeof(rz1)
    @test isone(irz1 * rz1)
    @test one(rz1) isa typeof(rz1)
    @test isa((irz1 * rz1), typeof(rz1))
end
