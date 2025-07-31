@testset "matrices" begin
    m1 = Matrix2x2(0.6233153827930635, 0.3551957329357007, 0.7085367960647345, 0.9084760155714717)
    m2 = Matrix2x2(0.04167333132967199, 0.6027145396563633, 0.42901908897608376, 0.7682825359251977)

    m1a = collect(m1)
    m2a = collect(m2)

    @test m1 == m1a
    @test m2 == m2a
    @test isapprox(m1 * m2, m1a * m2a)
    @test m1 + m2 == m1a + m2a
end

@testset "gates" begin
    T = Int
    (x, y, z) = (Xgate(T), Ygate(T), Zgate(T))
    T = Domega{Int}
    (xd, yd, zd) = (Xgate(T), Ygate(T), Zgate(T))
    @test x * y == im * z
    @test y * x == -im * z
    @test y * z == im * x
    @test z * y == -im * x

    # @test xd * yd == im * zd
    # @test yd * xd == -im * zd
    # @test yd * zd == im * xd
    # @test zd * yd == -im * xd
end
