import LinearAlgebra.lu
include("./Laplace2D.jl")
export pasovnaMatrika, ZgornjePasovnaMatrika, SpodnjePasovnaMatrika, lu, polna, pasovnaLaplace, resiPasovna

"""
    pasovnaMatrika(k1::Int, k2::Int, diagonale::Array{T,2})

Podatkovni tip, ki predstavlja pasovno matriko: razpršena matrika, 
ki ima vse neničelne vrednosti poleg na glavni diagonali še na k1 spodnjih in k2 zgornjih diagonalah.

# Argumenti
- `k1::Int`: število spodnjih diagonal
- `k2::Int`: število zgornjih diagonal
- `diagonale::Array{T,2}`: matrika (n x k1+k2+1 dimenzij) diagonal razvrščenih v stolpce: 
 najprej k1 spodnjih diagonal, nato glavna diagonala, ki ji sledi še k2 zgornjih diagonal

# Primer
```jldoctest
julia> A = [0 0 0 3 4 1 1; 0 0 6 7 8 1 0; 0 9 10 11 12 0 0; 2 13 14 15 0 0 0]
4×7 Array{Int64,2}:
 0   0   0   3   4  1  1
 0   0   6   7   8  1  0
 0   9  10  11  12  0  0
 2  13  14  15   0  0  0

julia> P = pasovnaMatrika(3,3,A)
pasovnaMatrika{Int64}(3, 3, [0 0 … 1 1; 0 0 … 1 0; 0 9 … 0 0; 2 13 … 0 0])

julia> polna(P)
4×4 Array{Float64,2}:
 3.0   4.0   1.0   1.0
 6.0   7.0   8.0   1.0
 9.0  10.0  11.0  12.0
 2.0  13.0  14.0  15.0
```
""" 
struct pasovnaMatrika{T}
    k1::Int
    k2::Int
    diagonale::Array{T,2} #lahko tudi list vektorjev, da bo brez ničelnih elementov. potem vcat dodamo 0 v konstruktorju
    
    function pasovnaMatrika(k1::Int, k2::Int, diagonale::Array{T,2}) where T
        m, n = size(diagonale)
        if (k1+k2+1 != n)
            error("Dimenzije matrike se ne ujemajo s podanim stevilom diagonal k1 in k2")
        else
            new{T}(k1,k2,diagonale)
        end
    end 
end

"""
    ZgornjePasovnaMatrika(k::Int, diagonale::Array{T,2})

Podatkovni tip, ki predstavlja zgornje pasovno matriko: razpršena matrika, 
ki ima vse neničelne vrednosti poleg na glavni diagonali še na k zgornjih diagonalah.

# Argumenti
- `k::Int`: število zgornjih diagonal
- `diagonale::Array{T,2}`: matrika (n x k+1 dimenzij) diagonal razvrščenih v stolpce: 
 najprej glavna diagonala, ki ji sledi še k zgornjih diagonal

# Primer
```jldoctest
julia> B = [1 2 3;1 2 3;1 2 0;1 0 0]
4×3 Array{Int64,2}:
 1  2  3
 1  2  3
 1  2  0
 1  0  0

julia> Z = ZgornjePasovnaMatrika(2,B)
ZgornjePasovnaMatrika{Int64}(2, [1 2 3; 1 2 3; 1 2 0; 1 0 0])

julia> polna(Z)
4×4 Array{Float64,2}:
 1.0  2.0  3.0  0.0
 0.0  1.0  2.0  3.0
 0.0  0.0  1.0  2.0
 0.0  0.0  0.0  1.0
```
""" 
struct ZgornjePasovnaMatrika{T}
    k::Int
    diagonale::Array{T,2} #lahko tudi list vektorjev, da bo brez ničelnih elementov. potem vcat dodamo 0 v konstruktorju
    
    function ZgornjePasovnaMatrika(k::Int, diagonale::Array{T,2}) where T
        m, n = size(diagonale)
        if (k+1 != n)
            error("Dimenzije matrike se ne ujemajo s podanim stevilom diagonal k")
        else
            new{T}(k,diagonale)
        end
    end 
end

"""
    SpodnjePasovnaMatrika(k::Int, diagonale::Array{T,2})

Podatkovni tip, ki predstavlja spodnje pasovno matriko: razpršena matrika, 
ki ima vse neničelne vrednosti poleg na glavni diagonali še na k spodnjih diagonalah.

# Argumenti
- `k::Int`: število spodnjih diagonal
- `diagonale::Array{T,2}`: matrika (n x k+1 dimenzij) diagonal razvrščenih v stolpce: 
 najprej k spodnjih diagonal, ki ji sledi še glavna diagonala

# Primer
```jldoctest
julia> C = [0 0 1;0 2 1;3 2 1;3 2 1]
4×3 Array{Int64,2}:
 0  0  1
 0  2  1
 3  2  1
 3  2  1

julia> S = SpodnjePasovnaMatrika(2,C)
SpodnjePasovnaMatrika{Int64}(2, [0 0 1; 0 2 1; 3 2 1; 3 2 1])

julia> polna(S)
4×4 Array{Float64,2}:
 1.0  0.0  0.0  0.0
 2.0  1.0  0.0  0.0
 3.0  2.0  1.0  0.0
 0.0  3.0  2.0  1.0
```
""" 
struct SpodnjePasovnaMatrika{T}
    k::Int
    diagonale::Array{T,2} #lahko tudi list vektorjev, da bo brez ničelnih elementov. potem vcat dodamo 0 v konstruktorju
    
    function SpodnjePasovnaMatrika(k::Int, diagonale::Array{T,2}) where T
        m, n = size(diagonale)
        if (k+1 != n)
            error("Dimenzije matrike se ne ujemajo s podanim stevilom diagonal k")
        else
            new{T}(k,diagonale)
        end
    end 
end

"""
Pomožni podatkovni tip, ki predstavlja unijo vseh tipov pasovnih matrik (pasovnaMatrika, ZgornjePasovnaMatrika, SpodnjePasovnaMatrika). 
"""
pasovne{T} = Union{pasovnaMatrika{T}, ZgornjePasovnaMatrika{T}, SpodnjePasovnaMatrika{T}}

"""
    e = getindex(P::pasovnaMatrika, i::Int, j::Int)

Vrne vrednost elementa pasovne matrike P[i, j].

# Argumenti
- `P::pasovnaMatrika`: pasovna matrika
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Rezultat
e::T je vrednost izbranega elementa

# Primer
```jldoctest
julia> polna(P)
4×4 Array{Float64,2}:
 3.0   4.0   1.0   1.0
 6.0   7.0   8.0   1.0
 9.0  10.0  11.0  12.0
 2.0  13.0  14.0  15.0

julia> P[3,2]
10
```
""" 
function Base.getindex(P::pasovnaMatrika, i::Int, j::Int)
    m, n = size(P.diagonale)
    d = P.k1 + 1 #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        return 0
    else
        return P.diagonale[x,y]
    end
end

"""
    e = getindex(P::ZgornjePasovnaMatrika, i::Int, j::Int)

Vrne vrednost elementa zgornje pasovne matrike P[i, j].

# Argumenti
- `P::ZgornjePasovnaMatrika`: zgornje pasovna matrika
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Rezultat
e::T je vrednost izbranega elementa

# Primer
```jldoctest
julia> polna(Z)
4×4 Array{Float64,2}:
 1.0  2.0  3.0  0.0
 0.0  1.0  2.0  3.0
 0.0  0.0  1.0  2.0
 0.0  0.0  0.0  1.0

julia> Z[3,2]
0
```
""" 
function Base.getindex(P::ZgornjePasovnaMatrika, i::Int, j::Int)
    m, n = size(P.diagonale)
    d = 1 #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        return 0
    else
        return P.diagonale[x,y]
    end
end

"""
    e = getindex(P::SpodnjePasovnaMatrika, i::Int, j::Int)

Vrne vrednost elementa spodnje pasovne matrike P[i, j].

# Argumenti
- `P::SpodnjePasovnaMatrika`: spodnje pasovna matrika
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Rezultat
e::T je vrednost izbranega elementa

# Primer
```jldoctest
julia> polna(S)
4×4 Array{Float64,2}:
 1.0  0.0  0.0  0.0
 2.0  1.0  0.0  0.0
 3.0  2.0  1.0  0.0
 0.0  3.0  2.0  1.0

julia> S[3,2]
2
```
""" 
function Base.getindex(P::SpodnjePasovnaMatrika, i::Int, j::Int)
    m, n = size(P.diagonale)
    d = n #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        return 0
    else
        return P.diagonale[x,y]
    end
end

"""
    setindex!(P::pasovnaMatrika, val::T, i::Int, j::Int)

Nastavi vrednost elementa pasovne matrike P[i, j].

# Argumenti
- `P::pasovnaMatrika`: pasovna matrika
- `val::T`: nova vrednost elementa
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Primer
```jldoctest
julia> P[1,1] = 100
100

julia> polna(P)
4×4 Array{Float64,2}:
 100.0   4.0   1.0   1.0
   6.0   7.0   8.0   1.0
   9.0  10.0  11.0  12.0
   2.0  13.0  14.0  15.0
```
""" 
function Base.setindex!(P::pasovnaMatrika, val::T, i::Int, j::Int) where T
    m, n = size(P.diagonale)
    d = P.k1 + 1 #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        error("Izbrani element ni diagonalni")
    else
        P.diagonale[x,y] = val
        return P
    end
end

"""
    setindex!(P::ZgornjePasovnaMatrika, val::T, i::Int, j::Int)

Nastavi vrednost elementa zgornje pasovne matrike P[i, j].

# Argumenti
- `P::ZgornjePasovnaMatrika`: zgornje pasovna matrika
- `val::T`: nova vrednost elementa
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Primer
```jldoctest
julia> Z[1,1] = 100
100

julia> polna(Z)
4×4 Array{Float64,2}:
 100.0  2.0  3.0  0.0
   0.0  1.0  2.0  3.0
   0.0  0.0  1.0  2.0
   0.0  0.0  0.0  1.0
```
""" 
function Base.setindex!(P::ZgornjePasovnaMatrika, val::T, i::Int, j::Int) where T
    m, n = size(P.diagonale)
    d = 1 #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        error("Izbrani element ni diagonalni")
    else
        P.diagonale[x,y] = val
        return P
    end
end

"""
    setindex!(P::SpodnjePasovnaMatrika, val::T, i::Int, j::Int)

Nastavi vrednost elementa spodnje pasovne matrike P[i, j].

# Argumenti
- `P::SpodnjePasovnaMatrika`: spodnje pasovna matrika
- `val::T`: nova vrednost elementa
- `i::Int`: indeks vrstice
- `j::Int`: indeks stolpca

# Primer
```jldoctest
julia> S[1,1] = 100
100

julia> polna(S)
4×4 Array{Float64,2}:
 100.0  0.0  0.0  0.0
   2.0  1.0  0.0  0.0
   3.0  2.0  1.0  0.0
   0.0  3.0  2.0  1.0
```
""" 
function Base.setindex!(P::SpodnjePasovnaMatrika, val::T, i::Int, j::Int) where T
    m, n = size(P.diagonale)
    d = n #indeks glavne diagonale
    y = d + (j - i) #d + indeks stranske diagonale
    x = i #element v diagonali

    if (y < 1 || y > n || x > m)
        error("Izbrani element ni diagonalni")
    else
        P.diagonale[x,y] = val
        return P
    end
end

"""
    firstindex(P::pasovne)
    
Vrne indeks prvega elementa pasovnih matrik P.

# Argumenti
P::pasovne: pasovna matrika

# Primer
```jldoctest
julia> firstindex(P)
1
```
""" 
Base.firstindex(P::pasovne) = 1

"""
    lastindex(P::pasovne)

Vrne indeks zadnjega elementa pasovnih matrik P.

# Argumenti
P::pasovne: pasovna matrika

# Primer
```jldoctest
#Polna P je 16x16 matrika
julia> lastindex(P)
16
""" 
Base.lastindex(P::pasovne) = size(P.diagonale,1)^2

"""
    r = *(P::pasovnaMatrika, v)

Množenje pasovne matrike P z desne z vektorjem v.

# Argumenti
- `P::pasovnaMatrika`: pasovna matrika dimenzij n x n
- `v::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
r::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia> P = pasovnaMatrika(3,3,[0 0 0 3 4 1 1; 0 0 6 7 8 1 0; 0 9 10 11 12 0 0; 2 13 14 15 0 0 0])
pasovnaMatrika{Int64}(3, 3, [0 0 … 1 1; 0 0 … 1 0; 0 9 … 0 0; 2 13 … 0 0])

julia> P*[1,2,3,4]
4×1 Array{Float64,2}:
 18.0
 48.0
 110.0
 130.0
```
""" 
function Base.:*(P::pasovnaMatrika, v)
    m, n = size(P.diagonale)
    r = zeros(m)

    for i=1:m
        for j=max(1,i-P.k1):min(i+P.k2,m)
            r[i] = r[i] + v[j]*P[i,j]
        end
    end

    return r 
end

"""
    r = *(P::ZgornjePasovnaMatrika, v)

Množenje zgornje pasovne matrike P z desne z vektorjem v.

# Argumenti
- `P::ZgornjePasovnaMatrika`: zgornje pasovna matrika dimenzij n x n
- `v::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
r::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia> Z = ZgornjePasovnaMatrika(2,[1 2 3;1 2 3;1 2 0;1 0 0])
ZgornjePasovnaMatrika{Int64}(2, [1 2 3; 1 2 3; 1 2 0; 1 0 0])

julia> Z*[1,2,3,4]
4×1 Array{Float64,2}:
 14.0
 20.0
 11.0
 4.0
```
""" 
function Base.:*(P::ZgornjePasovnaMatrika, v)
    m, n = size(P.diagonale)
    r = zeros(m)

    for i=1:m
        for j=i:min(i+P.k,m)
            r[i] = r[i] + v[j]*P[i,j]
        end
    end

    return r 
end

"""
    r = *(P::SpodnjePasovnaMatrika, v)

Množenje spodnje pasovne matrike P z desne z vektorjem v.

# Argumenti
- `P::SpodnjePasovnaMatrika`: spodnje pasovna matrika dimenzij n x n
- `v::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
r::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia> S = SpodnjePasovnaMatrika(2,[0 0 1;0 2 1;3 2 1;3 2 1])
SpodnjePasovnaMatrika{Int64}(2, [0 0 1; 0 2 1; 3 2 1; 3 2 1])

julia> S*[1,2,3,4]
4×1 Array{Float64,2}:
 1.0
 4.0
 10.0
 16.0
```
"""
function Base.:*(P::SpodnjePasovnaMatrika, v)
    m, n = size(P.diagonale)
    r = zeros(m)

    for i=1:m
        for j=max(1,i-P.k):i
            r[i] = r[i] + v[j]*P[i,j]
        end
    end

    return r 
end

"""
    x = \\(P::ZgornjePasovnaMatrika, b)

Resi enacbo Px = b z obratnim vstavljanjem.

# Argumenti
- `P::ZgornjePasovnaMatrika`: zgornje pasovna matrika dimenzij n x n
- `b::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
x::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia>  Z = ZgornjePasovnaMatrika(1,[100.0 2.0; 100.0 2.0; 100.0 2.0; 100.0 0.0])
ZgornjePasovnaMatrika{Float64}(1, [100.0 2.0; 100.0 2.0; 100.0 2.0; 100.0 0.0])

julia> b=Z*[1,2,3,4]
4×1 Array{Float64,2}:
 104.0
 206.0
 308.0
 400.0

julia> Z\\b
4×1 Array{Float64,2}:
 1.0
 2.0
 2.0
 4.0
```
""" 
function Base.:\(P::ZgornjePasovnaMatrika, b)
    m, n = size(P.diagonale)
    x = zeros(m)

    #obratno vstavljanje
    x[m] = b[m]/P[m,m];
    for i=m-1:-1:1
        for j=i+1:min(i+P.k,m)
            x[i] = x[i] + P[i,j]*x[j];
        end
        x[i] = (b[i]-x[i])/P[i,i];
    end

   return x
end

"""
    x = \\(P::SpodnjePasovnaMatrika, b)

Resi enacbo Px = b z direktnim vstavljanjem.

# Argumenti
- `P::SpodnjePasovnaMatrika`: spodnje pasovna matrika dimenzij n x n
- `b::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
x::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia> S = SpodnjePasovnaMatrika(1,[0.0 100.0;1.0 100.0;1.0 100.0;1.0 100.0])
SpodnjePasovnaMatrika{Float64}(1, [0.0 100.0; 1.0 100.0; 1.0 100.0; 1.0 100.0])

julia> b=S*[1,2,3,4]
4×1 Array{Float64,2}:
 100.0
 201.0
 302.0
 403.0

julia> S\\b
4×1 Array{Float64,2}:
 1.0
 2.0
 3.0
 4.0
```
""" 
function Base.:\(P::SpodnjePasovnaMatrika, b)
    m, n = size(P.diagonale)
    x = zeros(m)

    #direktno vstavljanje
    x[1] = b[1]/P[1,1];
    for i=2:m
        for j=max(1,i-P.k):i-1
            x[i] = x[i] + P[i,j]*x[j];
        end
    x[i] = (b[i]-x[i])/P[i,i];
   end

   return x
end

"""
    x = \\(P::pasovnaMatrika, b)

Reši enačbo Px = b z LU razcepom.

# Argumenti
- `P::pasovnaMatrika`: pasovna matrika dimenzij n x n
- `b::Array{T,1}`: vektor dimenzij n x 1

# Rezultat
x::Array{T,1}: vektor dimenzij n x 1

# Primer
```jldoctest
julia> P=pasovnaMatrika(1,1,[0.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 0.0])
pasovnaMatrika{Float64}(1, 1, [0.0 100.0 2.0; 1.0 100.0 2.0; 1.0 100.0 2.0; 1.0 100.0 0.0])

julia> b=P*[1,2,3,4]
4×1 Array{Float64,2}:
 104.0
 207.0
 310.0
 403.0

julia> P\\b
4×1 Array{Float64,2}:
 1.0
 2.0
 3.0
 4.0
```
""" 
function Base.:\(P::pasovnaMatrika, b)
    #uporabimo lu razcep
    L,U = lu(P)
    y = L\b
    x = U\y

    return x    
end

"""
    L, U = lu(P::pasovnaMatrika)

Naredi LU razcep pasovne matrike P.

# Argumenti
- `P::pasovnaMatrika`: pasovna matrika, ki jo želimo razcepiti

# Rezultat
L::SpodnjePasovnaMatrika: spodnje pasovna matrika razcepa
U::ZgornjePasovnaMatrika: zgornje pasovna matrika razcepa

# Primer
```jldoctest
julia> P=pasovnaMatrika(1,1,[0.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 2.0;1.0 100.0 0.0])
pasovnaMatrika{Float64}(1, 1, [0.0 100.0 2.0; 1.0 100.0 2.0; 1.0 100.0 2.0; 1.0 100.0 0.0])

L,U = lu(P)
(SpodnjePasovnaMatrika{Float64}(1, [0.0 1.0; 0.01 1.0; 0.010002 1.0; 0.010002 1.0]), ZgornjePasovnaMatrika{Float64}(1, [100.0 2.0; 99.98 2.0; 99.98 2.0; 99.98 0.0]))

julia> polna(L)
4×4 Array{Float64,2}:
 1.0   0.0       0.0       0.0
 0.01  1.0       0.0       0.0
 0.0   0.010002  1.0       0.0
 0.0   0.0       0.010002  1.0

julia> polna(U)
4×4 Array{Float64,2}:
 100.0   2.0    0.0    0.0
   0.0  99.98   2.0    0.0
   0.0   0.0   99.98   2.0
   0.0   0.0    0.0   99.98
```
""" 
function lu(P::pasovnaMatrika)
    m, n = size(P.diagonale)

    #preveri diagonalno dominantnost
    for i=1:m
        vsotaVrstice = 0
        for j=max(1,i-P.k1):min(i+P.k2,m)#1:n
            vsotaVrstice = vsotaVrstice + abs(P[i,j])
        end
        if abs(P[i,i]) < (vsotaVrstice-abs(P[i,i]))
            error("Matrika ni diagonalno dominantna") 
        end
    end

    #razcep
    A = pasovnaMatrika(P.k1, P.k2, copy(P.diagonale)) #kopiraj, da ne mutiras originalni P
    #L = SpodnjePasovnaMatrika(P.k1, hcat(zeros(m, P.k1),ones(m,1)))
    
    for k=1:m-1
        mi1 = min(k+P.k1,m)
        mi2 = min(k+P.k2,m)
        for i = k+1:mi1
            if (A[i,k] == 0)
                continue
            end
            A[i,k] = A[i,k]/A[k,k]
            for j = k+1:mi2
                if (A[k,j] == 0)
                    continue
                end  
                A[i,j] = A[i,j] - A[i,k]*A[k,j]
            end
        end
    end

    L = SpodnjePasovnaMatrika(A.k1, hcat(A.diagonale[:,1:A.k1],ones(m,1)))
    U = ZgornjePasovnaMatrika(A.k2, A.diagonale[:,A.k1+1:end])
    return L,U
end

"""
    M = polna(P)

Pomožna funkcija, ki vrne polno matriko M, ki vsebuje tudi ničelne vrednsoti pasovne matrike P.

# Argumenti
A::pasovnaMatrika: pasovna matrika

# Rezultat
M::Array{T,2}: polna matrika

# Primer
```jldoctest
julia> Z = ZgornjePasovnaMatrika(2,[1 2 3; 1 2 3; 1 2 0; 1 0 0])
ZgornjePasovnaMatrika{Int64}(2, [1 2 3; 1 2 3; 1 2 0; 1 0 0])

julia> polna(Z)
4×4 Array{Float64,2}:
 1.0  2.0  3.0  0.0
 0.0  1.0  2.0  3.0
 0.0  0.0  1.0  2.0
 0.0  0.0  0.0  1.0
```
""" 
function polna(P)
    m, n = size(P.diagonale)
    M = zeros(m, m)
    for i=1:m
        for j=1:m
            M[i,j] = getindex(P,i,j)
        end
    end

    return M
end

"""
    A = pasovnaLaplace(n, m)

Generira pasovno matriko A Laplaceovih enačb za reševanje minimalnih vpetih ploskev.

# Argumenti
- `n::Int`: prva dimenzija željene matrike
- `m::Int`: druga dimenzija željene matrike

# Rezultat
A::pasovnaMatrika: Laplaceova 2D matrika

# Primer
```jldoctest
julia> A = pasovnaLaplace(3, 3)
pasovnaMatrika{Float64}(3, 3, [0.0 0.0 … 0.0 1.0; 0.0 0.0 … 0.0 1.0; … ; 1.0 0.0 … 0.0 0.0; 1.0 0.0 … 0.0 0.0])

julia> polna(A)
9×9 Array{Float64,2}:
 -4.0   1.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0
  1.0  -4.0   1.0   0.0   1.0   0.0   0.0   0.0   0.0
  0.0   1.0  -4.0   0.0   0.0   1.0   0.0   0.0   0.0
  1.0   0.0   0.0  -4.0   1.0   0.0   1.0   0.0   0.0
  0.0   1.0   0.0   1.0  -4.0   1.0   0.0   1.0   0.0
  0.0   0.0   1.0   0.0   1.0  -4.0   0.0   0.0   1.0
  0.0   0.0   0.0   1.0   0.0   0.0  -4.0   1.0   0.0
  0.0   0.0   0.0   0.0   1.0   0.0   1.0  -4.0   1.0
  0.0   0.0   0.0   0.0   0.0   1.0   0.0   1.0  -4.0
```
""" 
function pasovnaLaplace(n, m)
    A = zeros(n*m,2*n+1)
    A[:,n+1] = ones(n*m,1)*-4 #glavna diagonala
    A[:,n] = ones(n*m,1) #prva spodnja
    for i=1:n:n*m
        A[i,n] = 0
    end
    A[:,n+2] = ones(n*m,1) #prva zgornja
    for i=n:n:n*m
        A[i,n+2] = 0
    end
    A[:,1] = vcat(zeros(n,1), ones(n*n-n,1)) #zadnja spodnja
    A[:,2*n+1] = vcat(ones(n*n-n,1), zeros(n,1)) #zadnja zgornja

    return pasovnaMatrika(n,n, A)
end

"""
    resiPasovna(::RobniProblemPravokotnik; nx=100, ny=100)

Izračunaj približek za rešitev robnega problema za operator 
z metodo deljenih diferenc. Uporabi pasovne matrike.

# Rezultat
- `Z::Matrix` je matrika vrednosti rešitve v notranjosti in na robu.
- `x::Vector` je vektor vrednosti abscise
- `y::Vector` je vektor vrednosti ordinate

# Primer

```julia
using Plots
robni_problem = RobniProblemPravokotnik(
    LaplaceovOperator(2),
    ((0, pi), (0, pi)),
    [sin, y->0, sin, y->0] 
)
Z, x, y = resiPasovna(robni_problem)
surface(x, y, Z)
```
"""
function resiPasovna(robni_problem; nx = 20, ny = 20)
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
    L = pasovnaLaplace(nx, ny) #uporabi pasovno matriko
    Z[2:end-1, 2:end-1] = reshape(L\b, nx, ny)
    return Z', x, y
end