# Overview

This folder contains the [Pironneau benchmark](https://doi.org/10.1017/S0022112074002023) for shape optimization in Stokes flow.

The code run in the article can be run by calling `make all`.

The file **Elasticity_solver.py** contains the mesh deformation scheme used in the second submission of the paper. It is a Linear-elasticity solver with pure
Neumann conditions.

The file **SingleMeshSolver.py** contains a single-mesh version of the Pironneau problem, using a linear elaticity equation for updating the mesh, and a steepest descent algorithm for optimization.

**StokesSolver.py** contains the multimesh implementation of the stokes
shape optimization problem.

**Stokes_Optimization.py** contains the version used in the first submission of
the paper, where we used the very complex [two stage eikonal deformation](https://arxiv.org/abs/1411.7663) scheme.