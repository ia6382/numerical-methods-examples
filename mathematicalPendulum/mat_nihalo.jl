import LinearAlgebra

#gravitacijska konstanta
const g = 9.80665

"""
dY = F(l, Y)

Vrni odvod funkcij.
Enačbi odvoda sta izpeljani iz diferencialne enačbe matematičnega nihala:
```math
    y_2^{\\prime} = \\frac{g}{l}sin(y_1)\\
    y_1^{\\prime} = y_2
```

# Parametra
l::Int64 dolžina vrvi nihala
Y::Array{Float64,1} vektor funkcijskih vrednosti [theta, dtheta]

# Rezultat
dY::Array{Float64,1} vektor odvoda funkcij  [dtheta, ddtheta]

# Primer
```jldoctest ;
julia> F(1, [1, 0.1])
2-element Array{Float64,1}:
  0.1
 -8.252011433166357
```
"""
function F(l, Y)
    dy1 = Y[2];
    dy2 = -(g/l)*sin(Y[1]);
    dY = [dy1, dy2];
    return dY;
end

"""
newY = rk_korak(h, l, prevY)

Ena iteracija/korak Runge-Kutta metode četrte stopnje.

# Parametri
l::Int64 dolžina vrvi nihala
h::Float64 širina podintervala
prevY::Array{Float64,1} vektor funkcijsih vrednosti v prejšnem koraku  [theta, dtheta]

# Rezultat
newY::Array{Float64,1} vektor funkcijsih vrednosti v novem koraku  [theta, dtheta]

# Primer
```jldoctest ;
julia> rk_korak(0.01, 1, [1, 0.1])
2-element Array{Float64,1}:
 1.0005873293604328
 0.01746068834694485
```
"""
function rk_korak(h, l, prevY)
    k1 = h*F(l, prevY);
    k2 = h*F(l, prevY+k1/2);
    k3 = h*F(l, prevY+k2/2);
    k4 = h*F(l, prevY+k3);
    newY = prevY + (k1+k2+k2+k3+k3+k4)/6;
    return newY;
end

"""
theta, Y, T = nihalo(l, t, theta0, dtheta0, n)

Izračunaj odmik matematičnega nihala v času t, pri podanih začetnih pogojih. Reši diferencialno enačbo:
```math
    {g\\over l}\\sin(\\theta(t))+\\theta^{\\prime\\prime}(t)=0, \\quad \\theta(0)=
    \\theta_0,\\ \\theta^{\\
    prime}(0)=\\theta^{\\prime}_0
```

# Parametri
l::Int64 dolžina vrvi nihala
t::Float64 čas v katerem želimo izračunati odmik nihala
theta0::Float64 začetni odmik nihala (v radianih)
dtheta0::Float64 začetna kotna hitrost nihala
n::Int64 število podintervalov

# Rezultat
theta0::Float64 odmik nihala ob času t (v radianih)
Y::Array{Float64,2} matrika funkcijskih vrednosti ob posameznih korakih  [theta_1 dtheta_1; ... ;theta_n dtheta_n]
T::Array{Float64,1} vektor časovnih korakov

# Primer
```jldoctest ;
julia> nihalo(1, 2, 1, 0.01, 2)
(-0.32085320653958393, [1.0 0.01; -0.320853 -1.25612], [0.0, 1.0])
```
"""
function nihalo(l, t, theta0, dtheta0, n)
    h = t/n;
    Y = zeros(n, 2);
    T = zeros(n);
 
    Y[1,:] = [theta0, dtheta0];
    T[1] = 0;

    for i=2:n
        Y[i, :] = rk_korak(h, l, Y[i-1, :])
        T[i] = T[i-1] + h;
    end
    theta = Y[n,1];

    return theta, Y, T;
end

"""
t = nihajni_cas(l, theta0, dtheta0, h, e)

Izračunaj nihajni čas s podanimi začetnimi pogoji. Za večjo natančnost na začetnem približku uporabi še Newtonovo metodo. 
Iščemo ničlo funkcije: 
```math
    \\theta (t) = \\theta_0
    \\theta (t) - \\theta_0 = 0
```

# Parametri
l::Int64 dolžina vrvi nihala
theta0::Float64 začetni odmik nihala (v radianih)
dtheta0::Float64 začetna kotna hitrost nihala
h::Float64 širina podintervala
e::Float64 meja za natančnost Newtonove metode

# Rezultat
t::Float64 nihajni čas (perioda)

# Primer
```jldoctest ;
julia> nihajni_cas(1, 0, 3.14, 0.01, 1e-10)
2.1541773466616854
```
"""
function nihajni_cas(l, theta0, dtheta0, h, e)
    times_to_cross = 2;
    if (dtheta0 != 0)
        times_to_cross = 2;
    end
    t = 0;
    Y = [theta0, dtheta0];
    prevY = Y;

    #računamo dokler ne prečkamo začetni odmik
    while (true)
        t = t + h;
        Y = rk_korak(h, l, Y)

        if ((prevY[1] > theta0 && Y[1] <= theta0) || (prevY[1] < theta0 && Y[1] >= theta0))
            if (times_to_cross == 1)
                break;
            else
                times_to_cross = times_to_cross -1;
            end
        end

        prevY = Y;
    end
    #println("približek ", t);

    #izboljšamo približek z Newtonovo metodo
    prevt = t - h;
    z = prevY;
    n = 1;
    while (abs(z[1]-theta0) > e && n < 100)
        t = prevt - (z[1]-theta0)/z[2];
        h = t - prevt;
        z = rk_korak(h, l, z);

        prevt = t;
        n = n+1;
    end

    return t;
end

"""
y = harm_nihalo(l, t, theta0, fi, n)

Aproksimacija matemtaičnega nihala, ki se za majhne odmike obnaša kot preprosto harmonično gibanje.

# Parametri
l::Int64 dolžina vrvi nihala
t::Float64 čas do katerega merimo kotne odmike nihala
theta0::Float64 začetni odmik nihala (v radianih)
fi::Float64 konstanta nihanja
n::Int64 število korakov

# Rezultat
Y::Array{Float64,1} vektor odmikov nihala ob posameznih korakih

# Primer
```jldoctest ;
julia> harm_nihalo(1, 3, 1, 0, 3)
3-element Array{Float64,1}:
  1.0
 -0.9999496444620966
  0.9997985829197464
```
"""
function harm_nihalo(l, t, theta0, fi, n)
    h = t/n;
    y = zeros(n);

    w = sqrt(g/l);
    t = 0;
    for i=1:n
        y[i] = theta0*cos(w*t + fi);
        t = t + h;
    end

    return y;
end


"""
Pomožna funkcija za analiziranje periode, ki sem jo uporabil za generiranje grafov dokumentacije. Ni del naloge. 
"""
function analiziraj_periodo(n)
    t = zeros(n);
    p = zeros(n);

    p[1] = 0.1;
    t[1] = nihajni_cas(1, 0, p[1], 0.01, 1e-10);
    for i=2:n
        p[i] = p[i-1] + 0.1
        t[i] = nihajni_cas(1, 0, p[i], 0.01, 1e-10); #glede na kotno hitrost
        #t[i] = nihajni_cas(1, p[i], 0.1, 0.01, 1e-10); #gleda na odmik
    end

    return t, p;
end

"""
theta0, Y, T = nihalo(1, 12, 3, 1, 300)
plot(Y[:, 1], Y[:, 2], label=[""], xlabel="odmik", ylabel="hitrost", title="theta = 3, dtheta = 1")
"""