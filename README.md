## Code used in MultiMesh Shape Opt article
Can be found at [arXiv](https://arxiv.org/pdf/1806.09821.pdf)

## Requirements
- Dolfin 2018.1.0
- Femorph (dokken/restructuring)
- IPOPT

## Installation
A docker container can be built with all the dependencies using the command
```
cd dockerfiles
docker build --tag mmshapeopt .
cd ..
```
To run the enviroment with current directory as shared with
```
docker run --name=mmshapeopt -ti -v $(pwd):/home/fenics/shared mmshapeopt
```
At an later instance, the container can be started with the following command
```
docker container start mmshapeopt
```

## Overview
The following examples can be found in the second submisson of the article:

- [Heat minimization of current carring multicables](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_MultiCable)
- [Drag minimization for obstacle in Stokes flow](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Stokes_Pironneau)
- [Drag minimization for orientation of nine obstacles in Stokes flow](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Stokes_rotation) 


The following examples were used in the first submission of the article:

- [Heat minimization problem for orientation of one obstacle](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_rotation)


The [Poisson comparasion](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_comparasion) folder is a folder with visual comparasion of the gradients for the shape derivatives using the Hadamard formulas for the MultiMesh FEM and traditional FEM.