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
docker run -ti -v $(pwd):/home/fenics/shared mmshapeopt
```