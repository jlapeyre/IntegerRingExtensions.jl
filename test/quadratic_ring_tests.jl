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

@testset "DRoot2" begin
    # x = ZRoot2(3)
    # @test x == 3
    # @test Int(x) == 3
    # r2 = root_two(ZRoot2{Int})
    # @test float(r2) == sqrt(2)
end
