using Test
using Random
using LinearAlgebra
import Base.*
using SparseArrays
include("precond_grad.jl")

@testset "Nepopolni Cholesky razcep" begin
    Random.seed!(0); D = rand(10, 10);
    L = tril(D);
    D = D'*D; 

    @test nep_chol(sparse(D)) ≈ cholesky(D).L
    @test test_cholesky(sparse(D)) == 16.007804381473342 #maksimalna napaka
end

@testset "Metoda kunjugiranih gradientov s predpogojevanjem" begin
    Random.seed!(0); S = create_sparseSPD(10)
    b = [1,2,3,4,5,6,7,8,9,10];
    L = nep_chol(S);
    x,i = conj_grad(S, b, L);

    @test S*x ≈ b
    @test x ≈ S\b
    @test test_grad(S, b, L) == 2.611244553918368e-13 #maksimalna napaka
end