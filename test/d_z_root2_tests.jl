@testset "Zroot2" begin
    x = Zroot2(3)
    @test x == 3
    @test Int(x) == 3
    r2 = root_two(Zroot2{Int})
    @test float(r2) == sqrt(2)
end

@testset "Droot2" begin
    # x = Zroot2(3)
    # @test x == 3
    # @test Int(x) == 3
    # r2 = root_two(Zroot2{Int})
    # @test float(r2) == sqrt(2)
end
