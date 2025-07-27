using CyclotomicRings: DOmega, ZOmega, CyclotomicRing
using QuadraticRings: QuadraticRing, ZRoot2
using RootOnes: RootOne
using SingletonNumbers: Two, InvTwo
using Dyadics: Dyadic

using Test

import CyclotomicRings: mul_root_two

include("test_aqua.jl")

@testset "CyclotomicRings.jl" begin
    @testset "cyclotomic" begin
        d = DOmega(1,2,3,4)
        z = ZOmega(1,2,3,4)
        for i in 0:7
            r = RootOne{8}(i)
            romega = DOmega{Int}(r)
            @test r * d == romega * d
            @test r * z == romega * z
        end
        @test isapprox(float(mul_root_two(d)) / float(d), sqrt(2))
        @test isapprox(float(mul_root_two(z)) / float(z), sqrt(2))

    end

    @testset "cyclotomic and rootone" begin
        @test CyclotomicRing{4}(0,1,0,0) == RootOne{8}(1)
        @test CyclotomicRing{4}(1,0,0,0) == RootOne{8}(0)
        @test CyclotomicRing{4}(0,0,1,0) == RootOne{8}(2)
        @test CyclotomicRing{4}(0,0,-1,0) == RootOne{8}(6)
        @test CyclotomicRing{3}(0,0,1) == RootOne{6}(2)
    end

    @testset "DOmega" begin
        cr = CyclotomicRing(1,2,3,4)
        @test isa(cr, CyclotomicRing{4, Int64})
        cromega = DOmega((1,2,3,4))
        @test cromega === DOmega(1,2,3,4)
        @test typeof(cromega) === DOmega{Int}
        @test cromega === DOmega{Int}(cr)
        z = DOmega(1,2,3,Dyadic(1,1))
        @test typeof(z) === DOmega{Int}
        @test typeof(DOmega(1,2,3,4)) === DOmega{Int}

        @test cromega == cromega
        @test DOmega(1,0,0,0) == one(cromega)
        @test isone(one(cr))
        @test isone(one(cromega))

        azf = float(z * conj(z))
        @test isreal(azf)
        @test real(azf) > 0
    end

    @testset "cyclotomic conversion" begin
        d = DOmega(1,2,3,4)
        @test isa(d, DOmega)
        @test !isa(d, ZOmega)

        z = ZOmega(d)
        @test isa(z, ZOmega)
        @test !isa(z, DOmega)

        dd = DOmega(z)
        @test !isa(dd, ZOmega)
        @test isa(dd, DOmega)

        @test real(z) isa QuadraticRing
        @test imag(z) isa QuadraticRing
        @test isapprox(float(z), float(real(z)) + im * float(imag(z)))

        z = ZOmega(1,2,3,4)
        zbig = ZOmega{BigInt}(z)
        zbigc = Complex(real(zbig) , imag(zbig))

        bigdiff = abs(big(zbig) - big(zbigc))
        floatdiff = abs(float(zbig) - float(zbigc))
        @test iszero(floatdiff)
        # This depends on precision(BigFloat). So it's fragile
        @test bigdiff < 1e-50

        zd = DOmega{Int64}((Dyadic{Int64, Int64}(3, 3), Dyadic{Int64, Int64}(1, 3), Dyadic{Int64, Int64}(3, 2), Dyadic{Int64, Int64}(0, 0)))
        zr = real(zd)
        zi = imag(zd)
        zr1 = QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(3, 3), Dyadic{Int64, Int64}(1, 4))
        zi1 = QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(3, 2), Dyadic{Int64, Int64}(1, 4))
        @test zr === zr1
        @test zi === zi1
        @test zr != zi1

        @test isa(abs2(zd), QuadraticRing)
        @test float(abs2(zd)) == abs2(float(zd))

        zd2 = DOmega{Int64}((Dyadic{Int64, Int64}(3, 3), Dyadic{Int64, Int64}(1, 3), Dyadic{Int64, Int64}(3, 2), Dyadic{Int64, Int64}(5, 5)))
        @test DOmega(complex(zd2)) === zd2
        # Lengths 5 and 7 are not supported well, and may not make sense as rings.
        # Might want to detect these and disallow
        # @test isa(CyclotomicRing{5, Dyadic{Int,Int}}((1,2,3//8,4,5)), CyclotomicRing{5, Dyadic{Int,Int}})
        # @test isa(CyclotomicRing{7, Dyadic}((1//16,2,3//8,4,5,6,7)), CyclotomicRing{7, Dyadic{Int,Int}})

    end

    @testset "cyclotomic D" begin
        z5 = CyclotomicRing{5}(0,1,0,0,0)
        @test float(z5) == float(RootOne{10}(1))

        z5 = CyclotomicRing{5}(0,0,1,0,0)
        @test float(z5) == float(RootOne{10}(2))

        z5 = CyclotomicRing{5}(0,0,0,1,0)
        @test float(z5) == float(RootOne{10}(3))

        z5 = CyclotomicRing{5}(0,0,0,0,1)
        @test float(z5) == float(RootOne{10}(4))
    end

    @testset "cyclotomic storage" begin
        z = DOmega(1,2,3,Dyadic(5,2))
        @test isbits(z)
        @test isa(z, DOmega{Int})
        zb = CyclotomicRing{4}(1,2,3,Dyadic(big(5),2))
        @test !isbits(zb)
        @test isa(zb, DOmega{BigInt})

        z2 = DOmega(Int128(1),2,3,Dyadic(5,2))
        @test isbits(z2)
        @test isa(z2, CyclotomicRing{4, Dyadic{Int128,Int}})
    end

    @testset "CyclotomicRing construction" begin
        @test ZOmega(3) isa ZOmega{Int}
        @test DOmega(10) isa DOmega{Int}
        @test DOmega{Int32}(10) isa DOmega{Int32}
        # These test both construction and conversion
        @test Int(ZOmega(3)) === 3
        @test Int(DOmega(3)) === 3

        @test typeof(DOmega(3//8)) == DOmega{Int}

        @test_throws MethodError DOmega(1,2,3)
        @test_throws MethodError ZOmega(1,2,3)

        @test ZRoot2(ZOmega(ZRoot2(1, 2))) === ZRoot2(1,2)
    end

    @testset "CyclotomicRing mixed artithemtic" begin
        x = DOmega(1,3//2, -5//4, 1//8)
        @test Two * (InvTwo * x) === (Two * InvTwo) * x
    end

end
