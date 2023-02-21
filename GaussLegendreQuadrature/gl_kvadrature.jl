import Base.*
using LinearAlgebra

"""
   T = gauss_quad(fun, a, b, n)

Izracunaj približek integrala funckije fun na intervalu [a, b].
Uporablja Gauss-Legendreovo integracijsko pravilo za interval [a, b]:
```math
\\int_{a}^{b}f(x)dx\\approx \\frac{b-a}{2} \\sum^n_{i=1}w_if(\\frac{b-a}{2}x_i+\\frac{a+b}{2})
```
Za dve točki Gauss-Legendrove kvardature sta uteži  w1 = w2 = 1 in vozlišči x1 = -1/sqrt(3) in x2 = -1/sqrt(3).
Po sestavljenem pravilu interval razdeli na n podintervalov in njihove približke nato sešteje v skkupni približek T.

# Parametri
fun funkcija, ki jo želimo integrirati
a::Int64 spodnja meja intervala
b::Int64 zgornja meja intervala
n::Int64 število podintervalov

# Rezultat
T::Float64 približek integrala funkcije fun

# Primer
```jldoctest ;
julia> gauss_quad(x->7*(x^3), 0, 5, 1)
1093.7500000000002
```
"""
function gauss_quad(fun, a, b, n)
    h = (b-a)/n;
    x = [-1/sqrt(3), 1/sqrt(3)];
    w = [1, 1];

    T = 0;
    for i=1:n
        #meje podintevala
        na = a+h*(i-1);
        nb = a+h*i;
        #enacba za Gauss-Legendreova dveh tock na intevalu [na, nb]
        I = w[1]*fun(((nb-na)/2)*x[1]+((na+nb)/2)) + w[2]*fun(((nb-na)/2)*x[2]+((na+nb)/2));
        I = I*((nb-na)/2);
        #sestavljeno pravilo
        T = T + I;
    end

    return T;
end

"""
    N, E = test_gq(fun, f, a, b)

Izracunaj napako Gauss-Legendreove integracije funkcije fun na intervalu [a, b] z različnim številom podintervalov.

# Parametri
fun funkcija, ki jo želimo integrirati
f::Float64 pravilna vrednost integrala
a::Int64 spodnja meja intervala
b::Int64 zgornja meja intervala

# Rezultat
N::Array{Int64,1} vektor števil podintervalov
E::Array{Float64,1} vektor napak integracije za določene podintervale

# Primer
```jldoctest ;
julia> f = 8495/12;

julia> fun(x) = 7*(x^3) - 8*(x^2) - 3*x - 3;

julia> julia> N,E =  test_gq(fun, f, 0, 5)

([1.0, 10.0, 100.0, 1000.0], [1.13687e-13, 1.13687e-13, 1.13687e-13, 2.27374e-13])

julia> scatter(N, E, xlabel="stevilo podintervalov", ylabel="napaka", scale=:log10)
```
"""
function test_gq(fun, f, a, b)
    N = zeros(8);
    E = zeros(8);
    
    N[1] = 1;
    E[1] = abs(f - gauss_quad(fun, a, b, 1));
    for i=2:8
        N[i] = 10*N[i-1];
        I = gauss_quad(fun, a, b, N[i])
        E[i] = abs(f - I);
    end

    return N, E;
end

"""
testni funkciji na intervalu [0, 5]
f::Float64 točna vrednost integrala funkcije 
"""
#funkcija 1
#f1 = 1.549931244944674137274408400730639012183184893966372210477969710681487208951511074986007223927691325;
#fun1(x) = sin(x)/x;
#funkcija 2
#f2 = 8495/12;
#fun2(x) = 7*(x^3) - 8*(x^2) - 3*x - 3;