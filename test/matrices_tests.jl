using LinearAlgebra: eigvecs, norm, dot
using IsApprox: Approx, isunitary

@testset "matrices" begin
    m1 = Matrix2x2(0.6233153827930635, 0.3551957329357007, 0.7085367960647345, 0.9084760155714717)
    m2 = Matrix2x2(0.04167333132967199, 0.6027145396563633, 0.42901908897608376, 0.7682825359251977)

    m1a = collect(m1)
    m2a = collect(m2)

    @test m1 == m1a
    @test m2 == m2a
    @test isapprox(m1 * m2, m1a * m2a)
    @test m1 + m2 == m1a + m2a

    m = Matrix2x2(1,2,3,4)
    mm = Matrix(m)
    for p in (0, 1, 2, 3, 10, 11, 20, 21)
        @test m^p == mm^p
    end
end

@testset "some matrix integraton tests" begin
    mh = Matrix2x2(Gate1{:H})
    @test isunitary(mh, Approx())
    (v1, v2) = columns(eigvecs(mh))
    @test isone(norm(v1), Approx())
    @test isone(norm(v2), Approx())
    @test abs(dot(v1, v2)) < 1e-14
end

@testset "SU2" begin
    ub = rand(SU2B)
    uc = SU2C(ub)
    ua = SU2(ub)
    mb = Matrix2x2(ub)
    mc = Matrix2x2(uc)
    ma = Matrix2x2(ua)
    @test isapprox(mb, mc)
    @test isapprox(mb, ma)

    # Broken
#    uc = rand(SU2C)
    # ub = SU2B(uc)
    # ua = SU2B(uc)
    # mb = Matrix2x2(ub)
    # mc = Matrix2x2(uc)
    # ma = Matrix2x2(ua)
    # @test isapprox(mb, mc)
    # @test isapprox(mb, ma)
end

# FIXME. Y and Z are not implemented anymore
# @testset "gates" begin
#     T = Int
#     (x, y, z) = map(Gate1{gatename}(Matrix2x2{T})(), (:X, :Y, :Z))
#     (Xgate(T), Ygate(T), Zgate(T))
#     T = DOmega{Int}
#     (xd, yd, zd) = (Xgate(T), Ygate(T), Zgate(T))
#     @test x * y == im * z
#     @test y * x == -im * z
#     @test y * z == im * x
#     @test z * y == -im * x

#     @test xd * yd == im * zd
#     @test yd * xd == -im * zd
#     @test yd * zd == im * xd
#     @test zd * yd == -im * xd

#     m = x + y + z
#     md = xd + yd + zd

#     n = 31 + 17im

#     @test n * m == float(n * md)
# end
