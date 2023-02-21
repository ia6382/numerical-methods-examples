# Numerical Methods Examples

# About
Several numerical analysis algorithms implemented in the [Julia language](https://julialang.org/).
Homework for the Numerical Mathematics course at the University of Ljubljana, Faculty of Computer and Information Science. The course covered several important algorithms and approaches for numerically solving problems in the following topics: systems of linear equations, eigenvalues, systems of nonlinear equations, interpolation and approximation, integration, derivation, systems of differential equations. 
Homework focused on specific scenarios where we could use the acquired knowledge of these topics.

# Content
Each folder contains the source code of the algorithm, short automatic tests and extensive documentation in the Slovenian language. To briefly summarize the content:
* Solving a system of linear equations with LU decomposition optimised specifically for the diagonally dominant band matrices
* Solving a system of linear equations with the conjugate gradient method and the incomplete Cholesky factorization as its preconditioner
* Barycentric interpolation variant of Lagrange polynomial interpolation 
* Approximating the definite integral of a function with the Gauss–Legendre quadrature method 
* Defining and solving a system of differential equations with the Runge-Kutta method for the angle of the mathematical model of a pendulum

# Instructions
* start Julia REPL (in Visual Studio Code IDE with Julia extension the shortcut is alt+j)
* include("sourceFile.jl") to use functions of the source file
* include("test.jl") to run the automatic tests
