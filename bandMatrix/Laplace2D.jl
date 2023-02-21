"""
2. laboratorijske vaje
"""

import LinearAlgebra.diagm
import Base:Matrix
export matrika, desne_strani, resi, LaplaceovOperator, RobniProblemPravokotnik

"""
Podatkovni tip brez vrednosti, ki predstavlja matriko za laplaceov operator v 
d dimenzijah.

# Primer
    L = LaplaceovOperator{2}()

vrne vrednost tipa `LaplaceovOperator{2}`, ki predstavlja 
laplaceov operator v dveh dimenzijah. 
"""
# Napiši definicijo strukture LaplaceovOperator
struct LaplaceovOperator{d} end

# Napiš pomožno funkcijo LaplaceovOperator(d), ki vrne vrednost
# tipa LaplaceovOperator{d}
LaplaceovOperator(d) = LaplaceovOperator{d}()

"""
    L = matrika(n, m, LaplaceovOperator(2))

Zapiši matriko za Laplaceov operator v 2D na pravokotnem območju. 
Matrika `L` je matrika sistema enačb za diskretizirano laplaceovo enačbo

```math
u_{i-1,j}+u_{i,j-1} - 4u_{ij} + u_{i+1,j}+u_{i,j+1} = 0.
```
"""
function matrika(n, m, ::LaplaceovOperator{2})
    L = diagm(0=>-4*ones(n), 1=>ones(n-1), -1=>ones(n-1))
    I = diagm(0=>ones(n)) # identiteta
    A = zeros(n*m, n*m)
    for j=1:m
        k = ((j-1)*n+1):(j*n) # indeksi v j-tem bloku  
        A[k,k] = L
        if j < m
            A[k,k.+ n] = I
            A[k.+n, k] = I
        end
    end
    return A
end

"""
    desne_strani(s, d, z, l, LaplaceovOperator(2))

Izračunaj desne strani pri reševanju robnega problema za Laplaceovo enačbo v 2
dimenzijah.
# Argumenti
- `s::Vector`: robne vrednosti na spodnjem robu
- `d::Vector`: robne vrednosti na desnem robu
- `z::Vector`: robne vrednosti na zgornjem robu
- `l::Vector`: robne vrednosti na levem robu
""" 
function desne_strani(s, d, z, l, ::LaplaceovOperator{2})
    n = length(s)
    m = length(l)
    b = zeros(n*m)
    b[1:n] -= s # j = 1
    b[n:n:end] -= d # i = n
    b[end-n+1:end] -= z # j = m
    b[1:n:end-n+1] -= l # i = 1
    return b
end

# definiraj podatkovno strukturo za robni problem na kvadratu. 
# Podatkovna struktura naj vsebuje operator, meje pravokotnika in
# vektor funkcij za vsako stranico posebej 
"""
    RobniProblemPravokotnik(operator, ((a, b), (c, d)), [f_s, f_d, f_z, f_l])

Definiraj robni problem za enačbo z danim diferencialnim operatorjem
```math
\\mathcal{L} u(x,y) = 0
```
na pravokotniku ``[a, b]\\times[c, d]``, kjer so vrednosti na robu podane s 
funkcijami ``u(x, c) = f_s(x)``, ``u(b, y) = f_d(y)``, 
``u(x, d) = f_z(x)`` in ``u(a, y) = f_l(y)``. 
"""
struct RobniProblemPravokotnik
    operator # operator
    meje # oglišča pravokotnika
    rp # funkcije na robu 
end

"""
    resi(::RobniProglemPravokotnik; nx=100, ny=100)

Izračunaj približek za rešitev robnega problema za operator 
z metodo deljenih diferenc.

# Rezultat
- `Z::Matrix` je matrika vrednosti rešitve v notranjosti in na robu.
- `x::Vector` je vektor vrednosti abscise
- `y::Vector` je vektor vrednosti ordinate

# Primer

```julia
using Plots
robni_problem = RobniProblemPravokotnik(
    LaplaceovOperator{2},
    ((0, pi), (0, pi)),
    [sin, y->0, sin, y->0] 
)
Z, x, y = resi(robni_problem)
surface(x, y, Z)
```
"""
function resi(robni_problem; nx = 20, ny = 20)
    (a, b), (c, d) = robni_problem.meje
    Z = zeros(nx + 2, ny + 2)
    x = LinRange(a, b, nx + 2)
    y = LinRange(c, d, ny + 2)
    Z[:, 1] = robni_problem.rp[1].(x)
    Z[end, :] = robni_problem.rp[2].(y)
    Z[:, end] = robni_problem.rp[3].(x)
    Z[1, :] = robni_problem.rp[4].(y)
    b = desne_strani(Z[2:end-1, 1], Z[end, 2:end-1],
     Z[2:end-1, end], Z[1, 2:end-1], robni_problem.operator)
    L = matrika(nx, ny, robni_problem.operator)
    Z[2:end-1, 2:end-1] = reshape(L\b, nx, ny)
    return Z', x, y
end
