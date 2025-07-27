using QuadraticRings: conj_root_two, QuadraticRing, norm_root_two, isunit,
    ZRoot2, Dyadic, QuadraticRing2, DRoot2, root_two

using RootOnes: omega
using SingletonNumbers: RootTwo

using Test

@testset "QuadraticRing{2, Int}" begin
    # Don't follow this example in real code!
    # If D is not literal or `const`, performance degrades by orders of magnitude!
    D = 2
    q = QuadraticRing{D}(1, 1)
    @test q === QuadraticRing{D}(1, 1) # different ways to instantiate.
    @test iszero(zero(q))
    @test iszero(zero(QuadraticRing{D, Int}))
    @test iszero(q - q)
    @test q^2 == q * q
    @test q^5 == q * q * q * q * q
    @test norm_root_two(q^5) == norm_root_two(q)
    @test conj_root_two(q) * q == norm_root_two(q)

    q2 = QuadraticRing2(1, 2)
    @test norm_root_two(q2) == -7
    @test conj_root_two(q2) * q2 == norm_root_two(q2)
    @test conj_root_two(q2) * q2 == -7.0

    u1 = QuadraticRing{D}(1, 1)
    u2 = QuadraticRing{D}(1, 0)
    @test isunit(u1)  # Units are not unique
    @test isunit(u2)

    @test convert(Int, QuadraticRing{D}(3, 0)) == 3
    @test_throws ArgumentError convert(Int, QuadraticRing{D}(3, 1)) # Can't convert 3 + sqrt(2) to Int
    @test float(QuadraticRing{D}(3, 2)) == 3 + 2 * sqrt(2)
end

@testset "ZRoot2" begin
    x = ZRoot2(3)
    @test x == 3
    @test Int(x) == 3
    r2 = root_two(ZRoot2{Int})
    @test float(r2) == sqrt(2)

    cnt = 0
    n = 10^3
    rng = 10^7
    for i in 1:n
        a = rand(-rng:rng)
        b = rand(-rng:rng)
        x = QuadraticRing{2}(a, b)
        if sign(x) != sign(float(x))
            cnt += 1
        end
    end
    @test cnt == 0

    res = QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(1, 0), Dyadic{Int64, Int64}(1, 1)) - QuadraticRing{2, Dyadic{Int64, Int64}}(Dyadic{Int64, Int64}(0, 0), Dyadic{Int64, Int64}(1, 1))*im

    @test 1 - complex(omega^3) === res
end

@testset "Quadratic Construction" begin
    @test ZRoot2(RootTwo) == ZRoot2(0, 1)
    @test isa(ZRoot2{BigInt}(RootTwo), ZRoot2{BigInt})
end

# TODO: errors: imaginary(DRoot2)
@testset "DRoot2" begin
    x = DRoot2(1, Dyadic(3, 2))
    @test typeof(x) === DRoot2{Int64, Int64}
    @test typeof(big(x)) === QuadraticRing{2, BigFloat}
    @test typeof(big(big(x))) === BigFloat
end

@testset "QuadraticRings fixes" begin
    r = complex(omega)
    @test isa(r^3, Complex{<:QuadraticRing})
end

@testset "Dyadic" begin
    x = DRoot2(Dyadic(3,2), 5)
    @test RootTwo * x == DRoot2(10, Dyadic(3,2))
    @test x * RootTwo == DRoot2(10, Dyadic(3,2))
    @test x === InvRootTwo * (RootTwo * x)
    @test x === RootTwo * (InvRootTwo * x)
end

#@testset "DRoot2" begin
    # x = ZRoot2(3)
    # @test x == 3
    # @test Int(x) == 3
    # r2 = root_two(ZRoot2{Int})
    # @test float(r2) == sqrt(2)
#end
