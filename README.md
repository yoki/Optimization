Optimization
============

Fortran 90 non-linear equations solver with documentation using Brent's method and Powells modified Hybrid method. 


What it does
----------
* It solves a system of nonlinear equations (hybrid.f90).
* It solves a nonlinear equation with Brent's method (fzero.f90). The single equation problem should not be treated as one case of the system of equations problem. Brent's method is only applicable to the single equation model, but it won't be stuck in the local solution, unlike hybrid.f90. 
* It provides detailed documentation that explains the intuition behind the algorithm and every subroutine in the code. The documentation contains details about why your model doesn't converge and what to do. 
* As a bonus, it does not depend on any external library. 

Motivation
----------
Canned optimization packages are useful when it works. But when it doesn't, I can do very little. Error messages didn't make sense. I didn't have a clue for the reason of non-convergence. I realized that I have to open the can. 

The challenge of opening the can is not about code itself. It is about documentations which closely follows the code. Comments in the code are usually not enough because we need to understand the derivations of formulas and intuitions behind the formulas. This package provides serious documentation. Symbols and subroutine structure in the documentation exactly match with code. 


License
-------
MIT
