using Test
include("gl_kvadrature.jl")
using LinearAlgebra
import Base.*

@testset "Gauss-Legendrove kvadrature dveh tock" begin

    @testset "sin(x)/x, interval [0, 5]" begin
        f1 = 1.549931244944674137274408400730639012183184893966372210477969710681487208951511074986007223927691325;
        fun1(x) = sin(x)/x

        I = gauss_quad(fun1, 0, 5, 122);
        @test I ≈ f1
        e = abs(f1 - I)
        @test e <= 1e-10
    end

    @testset "polinom tretje stopnje, interval [0, 5]" begin
        f2 = 8495/12;
        fun2(x) = 7*(x^3) - 8*(x^2) - 3*x - 3

        I = gauss_quad(fun2, 0, 5, 1);
        @test I ≈ f2
        e = abs(f2 - I)
        @test e <= 1e-10
    end
end