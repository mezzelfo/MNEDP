# Finite Elements Method for PDEs on Parametrized Manifolds

In this project we develop a MATLAB code from scratch for solving linear elliptic, parabolic and hyperbolic equations on parametrized 2-manifolds using an adaptation of the classical FEM.

For example, we can solve the wave equation on a torus.

![Alt text](https://github.com/mezzelfo/MNEDP/blob/Manifolds/torus_wave.gif?raw=true "Title")

Technical details:
- The code, for now, uses only first-order polynomial Courant finite elements but it has been designed with higher order generalizations in mind.
- The parametrization chart is supposed rectangular and it's triangulated with a regular Delaunay triangulation.
- In the parabolic case, the time advancement is performed with the Crank-Nicolson method while in the hyperbolic one Newmark's method is used.