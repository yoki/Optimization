Optimization
============

Fortran 90 non-linear equations solver with documentation using Brent's method and Powells modified Hybrid method. 


What it does
----------
* It solves a system of nonlinear equations (hybrid.f90).
* It solves a nonlinear equation with Brent's method (fzero.f90).
* It provides detailed documentation that explains the intuition behind the algorithm and every subroutine in the code. The documentation contains details on why your model doesn't converge and what to do. 
* As a bonus, it does not depend on any external library. 

Motivation
----------
I spent enough time wondering why my model doesn't converge, using Matlab's fsolve or other black-box functions for optimization. It is one thing that my model doesn't converge, but it is quite another thing that I don't have a clue for the reason of non-convergence. Error messages did't make sense. I realized that I have to open the black-box. 

However, it wasn't as easy as I expected, as I couldn't connect the codes and textbooks on the computation. There are a lot of good books (e.g. Numerical Optimization by Jorge Nocedal and Stephen Wright) for optimization and non-linear equations solver. There are also some public domain or open source codes for optimization (e.g. MINPACK in Fortran or Apache Commons for Java). I want a package seriously documented with equations. In this package, symbols and subroutine structure in the documentation exactly match with code. So it is easier to read code. 


Licence
-------
MIT
