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

@testset "Angle" begin
    @test isapprox(Dar(1/4), pi/4)
    @test !isapprox(Dar(1/4), pi/5)
    @test isapprox(Dar(1//4), pi/4)
    @test float(Dar(1/4)) == pi/4
    @test float(Dar(3//5)) == 3*pi/5
    @test Dar(1//4) == Dar(1/4)
    @test Dar(1//4) == pi/4

    for ang in (1/5, 1//5, 3, 2.7)
        d = Dar(ang)
        for (f, fpi) in ((cos, cospi), (sin, sinpi),  (cis, cispi))
            @test f(d) == fpi(ang)
        end
    end

    mrand = () -> 20 * rand() - 10
    checkadd = function ()
        x = mrand()
        y = mrand()
        if  !isapprox(Dar(x + y),  (Dar(x) + Dar(y)))
            return true
        end
        return false
    end
    N = 10^6
    @test iszero(sum(checkadd() for _ in 1:N))

    checkrange = function ()
        d = Dar(mrand())
        (-1 <= d.x <= 1) && return false
        true
    end
    @test iszero(sum(checkrange() for _ in 1:N))
end
