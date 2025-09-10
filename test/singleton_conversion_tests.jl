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
    @test Rational(InvTwo) === 1//2
    @test Rational{Int}(InvTwo) === 1//2
    biginvtwo = Rational{BigInt}(InvTwo)
    @test biginvtwo == 1//2
    @test biginvtwo isa Rational{BigInt}
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

@testset "Singletons One" begin
    @test_throws MethodError One * "dog"
    @test_throws MethodError "zebra" * One
    @test_throws MethodError One^"zebra"
    @test_throws MethodError One^2.0

    for n in (0, 1, 2, -1, -2, 3, -3)
        @test One^n === One
    end

    @test One * One === One
    @test Zero * One === Zero
    @test One * Zero === Zero
    @test isone(One)
    @test ! iszero(One)
    @test One * 2.1 === 2.1
    @test 0.5 * One === 0.5
    for T in (Two, InvTwo, RootTwo, InvRootTwo)
        @test One * T === T
        @test T * One === T
    end
end
