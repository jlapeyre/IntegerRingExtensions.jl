@testset "Zroot2" begin
    x = Zroot2(3)
    @test x == 3
    @test Int(x) == 3
    r2 = root_two(Zroot2{Int})
    @test float(r2) == sqrt(2)

    cnt = 0
    n = 10^3
    rng = 10^7
    for i in 1:n
        a = rand(-rng:rng)
        b = rand(-rng:rng)
        x = QuadraticRing{2}(a, b)
        if cmpzero(x) != cmp(float(x), 0)
            cnt += 1
        end
    end
    @test cnt == 0
end

@testset "Droot2" begin
    # x = Zroot2(3)
    # @test x == 3
    # @test Int(x) == 3
    # r2 = root_two(Zroot2{Int})
    # @test float(r2) == sqrt(2)
end
