@testset "ZRot" begin
    for theta in (1/2, 1/4, 0.5678, Dar(1/4), Dar(1//4))

        rz = zrot(theta)
        @test isapprox(theta, get_theta(rz))
        ang = - 2 * angle(rz[1])
        @test isapprox(ang, get_theta(rz))

        # rz2 = zrotpi(theta/pi)
        # @test isapprox(theta, get_theta(rz2))
        # ang = - 2 * angle(rz2[1])
        # @test isapprox(ang, get_theta(rz2))

        # rz3 = zrothalfpi(theta/(2*pi))
        # @test isapprox(theta, get_theta(rz3))
        # ang = - 2 * angle(rz2[1])
        # @test isapprox(ang, get_theta(rz3))
    end
end
