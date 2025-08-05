@testset "RingMatrices" begin
    m = Domega(1, 2, 3, 4)
    for n in 0:7
        @test inv(RootOne{8}(n)) * (RootOne{8}(n) * m ) == m
    end
end
