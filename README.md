# SpritzFEM - v2
## FEM program for linear-static analysis

This project is part of the course of COMPUTATIONAL MECHANICS OF MATERIAL for the MSc in Mechanical Engineering.

Goal is to build an all-around Fem code dealing with linear-elasticity (and eventually non-linear analysis and damage simulation).

The main engine is built around isoparametric elements that are integrated via Gauss-Legendre quadrature, 1D and 2D solver are under development.
3D tethraedal elements have been implemented altough at the moment the code for pressur-like BCs isn't ready yet and therefore their usage is limited to the case where non homogenous Dirichlet BCs are enforced.

The code of the main engine is built on Python and resolve to other commercial open-source software for the meshing part and post-processing visualisation of results suchs as Gmsh and Paraview.
