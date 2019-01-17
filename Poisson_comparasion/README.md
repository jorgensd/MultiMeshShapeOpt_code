# Summary

This folder contains files for visual comparasion of the gradients for the
FEM and the MultiMesh FEM using the Hadamard formulas.
A reference solution (discretly consistent gradient) is computed with the
method of mappings/material derivatives.

# Run examples
To generate the meshes needed, run
`python3 create_multiple_meshes.py`
using meshio and the python interface for gmsh called pygmsh.
To run the traditional FEM example, run
`python3 SingleMeshPoisson.py`
and the MultiMesh FEM example
`python3 MultiMeshPoisson.py`.

# Output
Output of the shape gradients can be found in the output folder.


# The pyadjoint file
To use the pyadjoint capabilites to [automatically compute shape derivatives](https://github.com/jorgensd/MultiMeshShapeOpt_code/blob/master/Poisson_comparasion/SingleMeshPoisson_pyadjoint.py) for the traditional finite element method, the following [pyadjoint](https://bitbucket.org/dolfin-adjoint/pyadjoint/pull-requests/72) version has to be used. This is not supplied with the docker image, as this is work under progress for another article.