using Test
include("baricentricna_interpolacija.jl")
using LinearAlgebra
import Base.*

@testset "Baricentrična Lagrangeva interpolacija" begin

    @testset "1. funkcija. Interval [-1, 1]" begin
        fun1(x) = MathConstants.e^-x^2

        p = -0.4
        x = fun1(p);
        i = bar_inter(fun1, p, -1, 1, 50);
        e = abs(x-i)
        @test e <= 1e-6
    end

    @testset "2. funkcija. Interval [0, 10]" begin
        fun2(x) = sin(x)/x

        p = 3
        x = fun2(p);
        i = bar_inter(fun2, p, 0, 10, 131);
        e = abs(x-i)
        @test e <= 1e-6
    end
    
    @testset "3. funkcija. Interval [1, 3]" begin
        fun3(x) = x^2 - 2*x

        p = 1.5
        x = fun3(p);
        i = bar_inter(fun3, p, 1, 3, 110);
        e = abs(x-i)
        @test e <= 1e-6

        #oscilira
        i = bar_inter(fun3, p, 1, 3, 3);
        e = abs(x-i)
        @test e <= 1e-6
    end
end
