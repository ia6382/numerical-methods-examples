using Test
include("pasovnaMatrika.jl")

@testset "Pasovne matrike operacije" begin
    A = [0 0 0 3 4 1 1; 0 0 6 7 8 1 0; 0 9 10 11 12 0 0; 2 13 14 15 0 0 0]
    P = pasovnaMatrika(3,3,A)
    @test polna(P) == [3 4 1 1;6 7 8 1;9 10 11 12;2 13 14 15]

    B = [1 2 3;1 2 3;1 2 0;1 0 0]
    Z = ZgornjePasovnaMatrika(2,B)
    @test polna(Z) == [1 2 3 0;0 1 2 3;0 0 1 2;0 0 0 1]

    C = [0 0 1;0 2 1;3 2 1;3 2 1]
    S = SpodnjePasovnaMatrika(2,C)
    @test polna(S) == [1 0 0 0;2 1 0 0; 3 2 1 0; 0 3 2 1]

    @test firstindex(P) == 1
    @test firstindex(Z) == 1
    @test firstindex(S) == 1

    @test lastindex(P) == 16
    @test lastindex(Z) == 16
    @test lastindex(S) == 16

    x = [1,2,3,4]
    @test P*x ==  [18.0,48.0,110.0,130.0]
    @test Z*x ==  [14.0,20.0,11.0,4.0] 
    @test S*x ==  [1.0,4.0,10.0,16.0]

    #posredno testira tudi lu razcep, indeksiranje in "deljenja z leve" (\) za vse tipe pasovnih matrik
    P=pasovnaMatrika(1,1,[0.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 0.0])
    b=P*x
    @test P\b == x
end

@testset "Reševanje Laplaceove 2D enačbe s pasovno matriko" begin
robni_problem = RobniProblemPravokotnik(
    LaplaceovOperator(2),
    ((0, pi), (0, pi)),
    [sin, y->0, sin, y->0] 
)
Zp, xp, yp = resiPasovna(robni_problem) #s pasovno matriko
Z, x, y = resi(robni_problem) #s polno matriko

@test Zp ≈ Z
@test xp ≈ x
@test yp ≈ y
end
