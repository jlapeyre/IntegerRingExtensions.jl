@testset "cyclotomic" begin
    d = Domega(1,2,3,4)
    z = Zomega(1,2,3,4)
    for i in 0:7
        r = RootOne8(i)
        romega = Domega{Int}(r)
        @test r * d == romega * d
        @test r * z == romega * z
    end
end
