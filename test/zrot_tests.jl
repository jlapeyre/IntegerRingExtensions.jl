@testset "ZRot" begin
    for theta in (1/2, 1/4, 0.5678, Dar(1/4), Dar(1//4))
        rz = zrot(theta)
        @test isapprox(theta, get_theta(rz))
        ang = - 2 * angle(rz[1])
        @test isapprox(ang, get_theta(rz))
    end

    rz1 = zrot(Dar(1//8))
    rz2 = zrot(Dar(3//8))
    rz3 = rz1 * rz2
    @test get_theta(rz3) === Dar(1//2)
end
