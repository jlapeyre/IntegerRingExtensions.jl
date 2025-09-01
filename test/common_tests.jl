@testset "common" begin
    @test isunit(1)
    @test isunit(-1)
    @test !isunit(2)
    @test !isunit(3)

    @test isunit(1.0 + 2.0 * im)
    @test isunit(1.0 + 0.0 * im)
    @test isunit(-1.23456 + 2.0 * im)

    @test isunit(1 + 0im)
    @test isunit(-1 + 0im)
    @test isunit(0 + 1im)
    @test isunit(0 + -1im)
    @test !isunit(Complex(1,1))

    @test isunit(7//13)
    @test isunit(Complex(7//13, 5//17))
end
