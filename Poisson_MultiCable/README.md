# Minimization of heat in a multicable
This folder contains source code for solving the heat minimization problem of
finding the placement of N internal cables in a MultiCable.

The subfolder [single_mesh])https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_MultiCable/single_mesh) contains the traditional FEM implementation used for comparasion.

# Examples in article #

- equilateral.py (Includes timings of MultiMesh operations)
- five_cables.py
- plot_taylor.py

# Run examples
Use the commands in the Makefile to get the results used in the article,
i.e. `make equilateral` runs the equilateral multicable example.

