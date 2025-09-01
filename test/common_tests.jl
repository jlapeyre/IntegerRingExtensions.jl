@testset "common" begin
    @test isunit(1)
    @test isunit(-1)
    @test !isunit(2)
    @test !isunit(3)
end
