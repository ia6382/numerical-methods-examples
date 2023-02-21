import Base.*
using LinearAlgebra
"""
   x = chebyeshev_nodes(n)

Vrni n Čebiševih točk na intervalu [-1,1], izračunanih po enačbi:
```math
cos(\\frac{(2i-1)\\pi}{2n}); i=1, ...,n
```

# Parameter
n::Int64 željeno število Čebišejevih točk

# Rezultat
x::Array{Int64,1} vektor Čebišejevih točk

# Primer
```jldoctest ;
julia> chebyeshev_nodes(3)
3-element Array{Float64,1}:
  0.8660254037844387
  6.123233995736766e-17
 -0.8660254037844387
```
"""
function chebyeshev_nodes(n)
    x = zeros(n);
    for i=1:n
        #x[i] = cos((i*pi)/n);
        x[i] = cos(((2*i-1)*pi)/2n); #formula iz wikipedie
    end
    return x;
end

"""
   x = weights(n)

Vrne n uteži za baricentrično Lagrangevo interpolacijo.:
```math
\\lambda=\\begin{cases}
(-1)^{i} & 0< i <1\\
(-1)^{i}\\frac{1}{2} & i=0,n.
\\end{cases}
```

# Parameter
n::Int64 željeno število uteži

# Rezultat
x::Array{Int64,1} vektor uteži

# Primer
```jldoctest ;
julia> weights(3)
3-element Array{Float64,1}:
  0.5
 -1.0
  0.5
```
"""
function weights(n)
    x = zeros(n);
    x[1] = 0.5;
    for i=1:n-2
        x[i+1] = (-1)^i;
    end
    x[n] = (-1)^(n-1) * 0.5;
    return x;
end

"""
   L = bar_inter(fun, p, a, b, n)

Interpoliraj vrednost funkcije v točki p, na intervalu [a, b]. 
Uporablja baricentrično Lagrangevo interpolacijo:
```math
L(x)=\\begin{cases}
\\frac{\\sum\\frac{f(x_{j})\\lambda_{j}}{x-x_{j}}}{\\sum\\frac{\\lambda_{j}}{x-x_{j}}} & x\\not=x_{j}\\
f(x_{j}) & \\text{sicer}
\\end{cases}
```
Točke (vozlišča) so predvidoma definirane na intervalu [-1, 1] zato jih najprej preslikamo na željeni interval [a, b].

# Parametri
fun funkcija, ki jo interpoliramo
p::Float64 točka za katero želimo dobiti vrednost funkcija
a::Int64 spodnja meja intervala
b::Int64 zgornja meja intervala
n::stopnja polinoma s katerim interpoliramo funkcijo

# Rezultat
L::Float64 vrednost funkcije v točki p

# Primer
```jldoctest ;
julia> fun(x) = sin(x)/x;

julia> bar_inter(fun, 3, 0, 10, 10)
  0.0488699778806767
```
"""
function bar_inter(fun, p, a, b, n)
    x = chebyeshev_nodes(n);
    if (a != -1 && b != 1)
        x = map(x->lin_map(x, a, b), x)
    end
    y = map(fun, x);
    w = weights(n);

    #preveri ali je zeljena vrednost p, ze podana tocka
    t = findall(x->x==p, x)
    if size(t,1) > 0
        println("tocka ze podana");
        L = y[t[1]];
        return L;
    end

    #druga baricentricna formula
    j = w ./ (p .- x);
    L = sum(j .* y)/sum(j);
    return L;
end

"""
    N, E = test_inter(fun, p, a, b)

Izracunaj napako Baricentrične Lagrangeve interpolacije funkcije fun na intervalu [a, b] z različnimi stopnjami polinomov.

# Parametri
fun funkcija, ki jo interpoliramo
p::Float64 točka za katero želimo dobiti vrednost funkcija
a::Int64 spodnja meja intervala
b::Int64 zgornja meja intervala

# Rezultat
N::Array{Int64,1} vektor stopnenj polinomov
E::Array{Float64,1} vektor napak interpolacije za določene stopnje polinomov

# Primer
```jldoctest ;
julia> fun(x) = MathConstants.e^-x^2
julia> N, E = test_inter(fun, -0.4, -1, 1)
([1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1.0e6, 1.0e7], [0.147856, 9.07686e-5, 1.71531e-7, 1.80315e-10, 1.70863e-13, 5.55112e-16, 2.22045e-16, 4.44089e-16])
julia> scatter(N, E, xlabel="stopnja polinom", ylabel="napaka", scale=:log10)
```
"""
function test_inter(fun, p, a, b)
    x = fun(p);
    N = zeros(8);
    E = zeros(8);
    
    N[1] = 1;
    E[1] = abs(x - bar_inter(fun, p, a, b, Int(N[1])));
    for i=2:8
        N[i] = 10*N[i-1];
        I = bar_inter(fun, p, a, b, Int(N[i]));
        E[i] = abs(x-I);
    end

    return N, E;
end

"""
    lx = lin_map(x, c, d; a=-1, b=1)

Linearno preslikaj vrednost x iz intervala [a, b] na interval [c, d].

# Parametri
fun funkcija, ki jo interpoliramo
x::Float64 vrednost za katero želimo preslikati
a::Int64 spodnja meja izvornega intervala
b::Int64 zgornja meja izvornega intervala
c::Int64 spodnja meja ciljnega intervala
d::Int64 zgornja meja ciljnega intervala

# Rezultat
lx::Float64 preslikana vrednost x

# Primer
```jldoctest ;
julia> lin_map(-0.4, 0, 10, -1, 1)
3.0
```
"""
function lin_map(x, c, d, a=-1, b=1)
    #za splosen interval:
    lx = ((d-c)/(b-a))*(x-a)+c;
    return lx;
end

"""
testne funkcije
"""
#1. funkcija. Interval [-1, 1]
#fun1(x) = MathConstants.e^-x^2;
#2. funkcija. Interval [0, 10]
#fun2(x) = sin(x)/x;
#3. funkcija. Interval [1, 3]
#fun3(x) = x^2 - 2*x;