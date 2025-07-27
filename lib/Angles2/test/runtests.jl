using Angles2
using Test

@testset "Angles2" begin
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

    @test Dar(3) / pi === 3
    @test (Dar(1) // 3).c === 1//3
    @test Dar(1) == pi
    @test pi == Dar(1)
    @test Dar(1) + pi === Dar(2)

    @test zero(Dar{Float64}) === Dar(0.0)
    @test zero(Dar{Int}) === Dar(0)
    @test zero(Dar(3)) === Dar(0)

    @test one(Dar(2)) === 1
    @test one(Dar(3.1)) === 1.0

    @test Dar(1) - Dar(2) === Dar(-1)
    @test Dar(1) - Dar(2.0) === Dar(-1.0)
    @test Dar(3.1) - Dar(2.3) === Dar(3.1 - 2.3)
end

include("test_aqua.jl")
