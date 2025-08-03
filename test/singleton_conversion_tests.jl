using IntegerExtensions
using Test

@testset "Singletons RootImag" begin
    ST = RootImag
    # Note AbstractFloat(3im) gives Inexact error
    @test ! canconvert(ST, AbstractFloat)
    @test ! canconvert(ST, Rational)
    @test canconvert(ST, Complex)
    for T in (Int, Int32, Int8, Int128, BigInt)
        @test ! canconvert(ST, T)
        @test ! canconvert(ST, Rational{T})
        @test canconvert(ST, Complex{T})
        @test ! canconvert(ST, Complex{Rational{T}})
        @test ! canconvert(ST, float(T))
        @test canconvert(ST, Complex{float(T)})
    end
end

@testset "Singletons Imag" begin
    ST = Imag
    # Note AbstractFloat(3im) gives Inexact error
    @test ! canconvert(ST, AbstractFloat)
    @test ! canconvert(ST, Rational)
    @test canconvert(ST, Complex)
    for T in (Int, Int32, Int8, Int128, BigInt)
        @test ! canconvert(ST, T)
        @test ! canconvert(ST, Rational{T})
        @test canconvert(ST, Complex{T})
        @test canconvert(ST, Complex{Rational{T}})
        @test ! canconvert(ST, float(T))
        @test canconvert(ST, Complex{float(T)})
    end
end

@testset "Singletons Two" begin
    ST = Two
    @test canconvert(ST, AbstractFloat)
    @test canconvert(ST, Rational)
    @test canconvert(ST, Complex)
    for T in (Int, Int32, Int8, Int128, BigInt)
        @test canconvert(ST, T)
        @test canconvert(ST, Rational{T})
        @test canconvert(ST, Complex{T})
        @test canconvert(ST, Rational{T})
        @test canconvert(ST, Complex{Rational{T}})
        @test canconvert(ST, float(T))
        @test canconvert(ST, Complex{float(T)})
    end
end

@testset "Singletons InvTwo" begin
    ST = InvTwo
    @test canconvert(ST, AbstractFloat)
    @test canconvert(ST, Rational)
    @test canconvert(ST, Complex)
    for T in (Int, Int32, Int8, Int128, BigInt)
        @test ! canconvert(ST, T)
        @test canconvert(ST, Rational{T})
        @test ! canconvert(ST, Complex{T})
        @test canconvert(ST, Complex{Rational{T}})
        @test canconvert(ST, float(T))
        @test canconvert(ST, Complex{float(T)})
    end
end

@testset "Singletons RootTwo" begin
    for RT in (RootTwo, InvRootTwo)
        @test canconvert(RT, AbstractFloat)
        @test ! canconvert(RT, Rational)
        @test canconvert(RT, Complex)
        for T in (Int, Int32, Int8, Int128, BigInt)
            @test ! canconvert(RT, T)
            @test ! canconvert(RT, Rational{T})
            @test ! canconvert(RT, Complex{T})
            @test ! canconvert(RT, Complex{Rational{T}})
            @test canconvert(RT, float(T))
            @test canconvert(RT, Complex{float(T)})
        end
    end
end
