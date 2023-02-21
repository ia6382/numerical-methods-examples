import LinearAlgebra.norm
import LinearAlgebra.cholesky 
import LinearAlgebra.tril
import Base.*
import LinearAlgebra.I
import LinearAlgebra.diagind

using SparseArrays

"""
    A = create_sparseSPD(n)

Generira razprseno simetricno pozitivno definitno (SPD) nakljucno matriko

# Parameter
n::Int dimenzija matrike

# Rezultat
A::SparseMatrixCSC razprsena SPD matrika

# Primer
```jldoctest
julia> S = create_sparseSPD(3);
julia> Matrix(S)
3×3 Array{Float64,2}:
 1.2063    0.0       0.772486
 0.0       1.20024   0.719139
 0.772486  0.719139  0.846356
```
"""
function create_sparseSPD(n)
    A = sprand(n, n, 0.4);
    r = rand(n,1);

    L = tril(A);
    A = L'+L; # zagotovi simetricost
    A[diagind(A)] = (r.+(n/7)) # zagotovi pozitivno definitnost, ker je diagonalno dominantna
    return A;
end


"""
    L = nep_chol(A)

Naredi nepopolni cholesky razcep za razprseno matriko A.
```math
    A = L'*L
```
Postopek je enak popolnemu Cholesky razcepu, le da ignoriramo nicelne elemente

# Parameter
AA::SparseMatrixCSC razprsena SPD matrika, ki jo zelimo razcepiti

# Rezultat
L::SparseMatrixCSC spodnje trikotna razprsena SPD matrika

# Primer
```jldoctest
julia> S = sparse([1,2,3,2,3],[1,2,2,3,3],[1.53542, 0.709665, 0.62429, 0.62429, 1.32583]);
julia> nep_chol(S)
3×3 SparseMatrixCSC{Float64,Int64} with 4 stored entries:
  [1, 1]  =  1.23912
  [2, 2]  =  0.842416
  [3, 2]  =  0.741071
  [3, 3]  =  0.881274
```
"""
function nep_chol(AA)
    A = copy(AA);
    n = size(A,1);
    I = rowvals(A); #return a vector of the row indices of A

    for k = 1:n
        A[k,k] = sqrt(A[k,k]);

        r = nzrange(A,k); #return the range of indices to the structural nonzero values of a sparse matrix column
        rows = I[r];
        nzK = rows[rows .> k]; #indeksi vrstic nenicelnih, podiagonalnih el. v stolpcu k
        
        #izracunaj r vektor
        A[nzK,k] = A[nzK,k]./A[k,k]; 

        for j = (k+1):n
            r = nzrange(A,j);
            rows = I[r];
            nzJ = rows[rows .>= j]; #indeksi vrstic nenicelnih, podiagonalnih(vkljucno z diagonalnimi) el. v stolpcu j
            
            #izracunaj diagonalno podmatriko R
            A[nzJ,j] = A[nzJ,j] - A[nzJ,k].*A[j,k];
        end
    end

    L = tril(A)
    return L;
end

"""
    e = test_cholesky(A)

Testira Cholesky razcep na A matriki in vrne najvecjo napako, za element matrike, ki se najbolj razlikuje od pravilnega razcepa.
Razcep testiramo, tako da primerjamo originalno matriko A in matriko sestavljeno iz rezultata razcepa L: 
```math
    ost = L'*L-A;
````
# Parameter
A::SparseMatrixCSC razprsena SPD matrika, ki jo zelimo razcepiti

# Rezultat
e::Float najvecja napaka razcepa

# Primer
```jldoctest 
julia> S = sparse([1,2,3,2,3],[1,2,2,3,3],[1.53542, 0.709665, 0.62429, 0.62429, 1.32583]);
julia> test_cholesky(S)
0.5491858892576076
```
"""
function test_cholesky(A)
    L = nep_chol(A);

    ost = L'*L-A;

    maxerr = maximum(abs.(ost));
    return maxerr;
end


"""
    x, i = conj_grad(A, b, L; N=length(b), e=1e-10)

Resi linearni sistem enacb Ax = b, z metodo konjugiranih gradientov s predpogojevanjem. 
Kot predpogojevalnik se pricakuje uporaba nepopolnega Cholesky razcepa.

# Parametri
A::SparseMatrixCSC razprsena SPD matrika, ki jo zelimo razcepiti
b::Array{Float,1} desna stran sistema enacb(rezultat)
L::SparseMatrixCSC spodnje trikotna razprsena matrika - rezultat Cholesky razcepa A matrike
N::Int maksimalno število iteracij
e::Float natančnost, ki določa ustavitveni pogoj metode

# Rezultat
x::Array{Float,1} rezultat sistema enacb
i::Int stevilo porabljenih iteracij za izracun rezultata pri dosezeni natancnosti

# Primer
```jldoctest 
julia> S = sparse([1,2,3,2,3],[1,2,2,3,3],[1.53542, 0.709665, 0.62429, 0.62429, 1.32583]);
julia> L = nep_chol(S);
julia> b = [1.0, 2.0, 3.0];
julia> x, i = conj_grad(S, b, L);
julia> x
3×1 Array{Float64,2}:
 0.6512875955764547
 1.4130089940827486
 1.5973937948938246
julia> i
3
```
"""
function conj_grad(A, b, L; N=length(b), e=1e-10)
    x = zeros(size(A,1), 1);

    r = b - A*x;
    z = L\(L'\r);
    p = z;

    rsold = r'*z;
    iter = 0;
    for i = 1:N
        iter = iter + 1;

        Ap = A*p;
        alfa = rsold/(p'*Ap);
        alfa = alfa[1];
        x = x + alfa*p;
        r = r - alfa*Ap;
        if norm(r) < e
            break
        end
        z = L\(L'\r);
        rsnew = r'*z;
        beta = rsnew/rsold;
        beta = beta[1];
        p = z + beta*p;

        rsold = rsnew;
    end

    return x, iter;
end

"""
    x, i = conj_grad_orig(A, b; N=length(b), e=1e-20)

Resi linearni sistem enacb Ax = b, z metodo konjugiranih gradientov brez predpogojevanja. 
Enaka funkciji conj_grad le, da ne vsebuje parametra L (oz. L je indentiteta). 
"""
conj_grad_orig(A, b; N=length(b), e=1e-10) = conj_grad(A, b, I; N=N, e=e)

"""
    e = test_grad(A, b, L)

Testira resevanje sistema enacb Ax = b, z metodo konjugiranih gradientov s predpogojevanjem. 
Vrne najvecjo napako (ostanek), po spodnji razliki:
```math
    e = A*x - b;
````
# Parametri
A::SparseMatrixCSC razprsena SPD matrika, ki jo zelimo razcepiti
b::Array{Float,1} desna stran sistema enacb(rezultat)
L::SparseMatrixCSC spodnje trikotna razprsena matrika - rezultat Cholesky razcepa A matrike

# Rezultat
e::Float najvecja napaka resitve sistema

# Primer
```jldoctest 
julia> S = sparse([1,2,3,2,3],[1,2,2,3,3],[1.53542, 0.709665, 0.62429, 0.62429, 1.32583]);
julia> L = nep_chol(S);
julia> b = [1,2,3];
julia> test_grad(S, b, L)
1.3322676295501878e-15
```
"""
function test_grad(A, b, L)
    x,i = conj_grad(A, b, L);

    r = A*x - b;
    maxerr = maximum(abs.(r));
    return maxerr;
end

"""
    range, data = compare_convergance(n, e)

Primerja konvergenco resevanje nakljucnega sistema enacb Ax = b, z metodo konjugiranih gradientov s predpogojevanjem in brez predpogojevanja. 
Za lazjo predstavo, lahko nato narisemo graf natacnosti resitve v odvisnosti od stevila iteracij, potrebnih za to natacnost (glej primer). 

# Parameter
n::Int dimenzija matrike
e::Float natančnost, ki določa ustavitveni pogoj metode

# Rezultat
range::UnitRange{Int64} število iteracij (x os grafa)
data::2-element Array{Array{Float64,2},1} vektor, ki vsebuje napake za metodo s in brez predpogojevanja [precondE, nonPrecondE] (y os grafa)

# Primer
```jldoctest
julia> range, data = compare_convergance(100, 1e-14)
julia> scatter(range, data, label=["predp.", "brez"], xlabel="iteracije", ylabel="napaka", scale=:log10);
```
"""
function compare_convergance(n, e)
    range = 10:100

    A = create_sparseSPD(n);
    b = rand(n,1);

    L = nep_chol(A);

    precondE = zeros(range[end]-range[1]+1,1)
    nonPrecondE = zeros(range[end]-range[1]+1,1)
 
    i = 1;
    for N=range
        x,iter = conj_grad(A, b, L, N=N, e=e);
        r = A*x - b;
        precondE[i] = maximum(abs.(r));

        x,iter = conj_grad_orig(A, b, N=N, e=e);
        r = A*x - b;
        nonPrecondE[i] = maximum(abs.(r));

        i = i+1;
    end

    return range, [precondE, nonPrecondE]
end