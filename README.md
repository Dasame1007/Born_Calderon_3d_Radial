This repository contains the Julia code used to produce the numerical results given in the article _Title_ [[arXiv]](https://link-url-here.org) [[Journal]](https://link-url-here.org) written by Fabricio Macià, Cristóbal J. Meroño and Daniel Sánchez-Mendoza.

The Julia module [Born_Calderon_3d_Radial.jl](src/Born_Calderon_3d_Radial.jl) has all necessary functions to compute the eigenvalues of the Dirichlet-to-Neumann map associated to a radial potential on the unit ball and construct from them the associated Born approximation.

The Jupyter notebook [Example.ipynb](src/Example.ipynb) presents an example of how to use the module.

__Dependencies:__
- ArbNumerics v1.5.0,
- SpecialFunctions, Jacobi, FFTW, Plots.
