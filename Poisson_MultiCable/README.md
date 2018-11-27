# Examples in article #
equilateral.py

# Installation #
For the examples to run, the following commands has to be run to obtain
dolfin 2018.1.0 with a compatible IPOPT solver. It also installs femorph 2.0,
where the only component needed is the VolumeNormal, a "CG1"-representation
of the FacetNormal in dolfin. To run the examples use python3.
```
docker run -ti -v $(pwd):/home/fenics/shared quay.io/dolfinadjoint/pyadjoint:2018.1.0
git clone https://bitbucket.org/Epoxid/femorph.git
cd femorph
git checkout dokken/restructuring
sudo pip3 install .
cd ../shared
```
