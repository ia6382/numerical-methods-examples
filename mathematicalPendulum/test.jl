using Test
include("mat_nihalo.jl")

@testset "matematicno nihalo" begin
    @testset "odvod dif. enacbe nihanja" begin
        dY = F(1, [1, 0.1])
        @test dY ≈ [0.1; -8.252011433166357]
    end

    @testset "Runge-Kutta korak 4. stopnje" begin
        newY = rk_korak(0.01, 1, [1, 0.1])
        @test newY ≈ [1.0005873293604328; 0.01746068834694485]
    end

    @testset "resitev dif. enacbe nihanja" begin
        t0, Z, T = nihalo(1, 2, 1, 0.01, 2)
        @test t0 ≈ -0.32085320653958393
        @test Z ≈ [1.0 0.01; -0.32085320653958393 -1.256117400309804]
        @test T ≈ [0.0, 1.0]
    end

    @testset "racunanje periode" begin
        t = nihajni_cas(1, 0, 3.14, 0.01, 1e-10)
        @test t ≈ 2.1541773466616854
    end

    @testset "harmonicno nihalo" begin
        Y = harm_nihalo(1, 3, 1, 0, 3)
        @test Y  ≈ [1.0; -0.9999496444620966; 0.9997985829197464]
    end
end