@testset "common" begin
    @test isunit(1)
    @test isunit(-1)
    @test !isunit(2)
    @test !isunit(3)

    @test isunit(1.0 + 2.0 * im)
    @test isunit(1.0 + 0.0 * im)
    @test isunit(-1.23456 + 2.0 * im)
end
