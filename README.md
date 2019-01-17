## Code used in MultiMesh Shape Opt article
Can be found at [arXiv](https://arxiv.org/pdf/1806.09821.pdf)

## Overview
The following examples can be found in the second submisson of the article:

- [Heat minimization of current carring multicables](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_MultiCable)
- [Drag minimization for obstacle in Stokes flow](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Stokes_Pironneau)
- [Drag minimization for orientation of nine obstacles in Stokes flow](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Stokes_rotation) 


The following examples were used in the first submission of the article:

- [Heat minimization problem for orientation of one obstacle](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_rotation)


The [Poisson comparasion](https://github.com/jorgensd/MultiMeshShapeOpt_code/tree/master/Poisson_comparasion) folder is a folder with visual comparasion of the gradients for the shape derivatives using the Hadamard formulas for the MultiMesh FEM and traditional FEM.

## Installation

### Docker
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

To run the attached docker notebook, run:
`
docker run --name mmshapeoptnb -w /home/fenics -v $(pwd):/home/fenics/shared -d -p 127.0.0.1:8888:8888 mmshapeopt 'jupyter-notebook --ip=0.0.0.0'
`


### Manual installation
If you do not want to use docker, you need the following packages:

- [Dolfin 2018.1.0](https://bitbucket.org/fenics-project/dolfin/src/2018.1.0.post2/) 
- Femorph, branch [dokken/restructuring](https://bitbucket.org/Epoxid/femorph/src/c7317791c8f00d70fe16d593344cb164a53cad9b/?at=dokken%2Frestructuring), which is pip installable
- IPOPT (v. 3.12.9)
- scipy (>=1.1.0)
- gmsh (v. 3.0.6)
- meshio, pygmsh
