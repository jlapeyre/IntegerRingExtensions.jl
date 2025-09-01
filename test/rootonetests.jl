@testset "RootOne" begin
    # Test that `angle` gives the same result obtained by converting to
    # float, and then extracting the angle.
    for i in 0:20
        root = RootOne{8}(i)
        @test RootOne{8}(i) === Omega(i)
        @test isapprox(angle(root), angle(float(root)))
        @test isone(root * inv(root))
        @test isone(inv(root) * root)
        @test isone(conj(root) * root)
        @test conj(root) == adjoint(root)
    end
    @test RootOne{8}() === RootOne{8}(1)
    r = RootOne{8}(1)
    @test isa(angle(r), Float64)
    @test isa(angle(BigFloat, r), BigFloat)

    # Test that BigFloat angle has full precision
    r5 = RootOne{8}(5)
    z = cis(angle(BigFloat, r5))
    @test real(z^2) < 1 / big(10)^75

    @test RootOne{8}(4) / RootOne{8}(2) === RootOne{8}(2)

    @test RootOne{8}(2) == RootOne{4}(1)
    @test RootOne{5}(0) == RootOne{17}(0)

    @test imaginary(RootOne{4}) == RootOne{4}(1)
    @test imaginary(RootOne{8}) == RootOne{8}(2)
    @test imaginary(RootOne{12}) == RootOne{12}(3)
    @test_throws ArgumentError imaginary(RootOne{6})

    @test sqrt_imaginary(RootOne{8}) == RootOne{8}(1)
    @test sqrt_imaginary(RootOne{16}) == RootOne{16}(2)
    @test_throws ArgumentError sqrt_imaginary(RootOne{12})

    @test CyclotomicRing{4}(0,1,0,0) == RootOne{8}(1)
    @test CyclotomicRing{4}(1,0,0,0) == RootOne{8}(0)
    @test CyclotomicRing{4}(0,0,1,0) == RootOne{8}(2)
    @test CyclotomicRing{4}(0,0,-1,0) == RootOne{8}(6)
    @test CyclotomicRing{3}(0,0,1) == RootOne{6}(2)

    for N in (1, 2, 4, 5, 8, 12, 16)
        for i in 0:(N-1)
            r = RootOne{N}(i)
            @test isreal(r) == isreal(float(r))
            @test isinteger(r) == isinteger(float(r))
        end
    end

    @test Integer(RootOne{2}(0)) === 1
    @test Integer(RootOne{2}(1)) === -1
    @test Int32(RootOne{2}(0)) === Int32(1)
    @test Int32(RootOne{2}(1)) === Int32(-1)
end
