# Test for verifying that T("-") is the top mesh value in dI
# Note: This tests only passes with following setup:
"""
git+https://bitbucket.org/magneano/ufl.git@merge-master
git+https://bitbucket.org/fenics-project/fiat.git@master
git+https://bitbucket.org/fenics-project/dijitso.git@master
git+https://bitbucket.org/magneano/ffc.git@evaluate-reference-basis
"""

from dolfin import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace
import numpy

def convergence_rates(E_values, eps_values):
    from numpy import log
    r = []
    for i in range(1, len(eps_values)):
        r.append(log(E_values[i]/E_values[i-1])/log(eps_values[i]/
                                                    eps_values[i-1]))

    print("Computed convergence rates: {}".format(r))
    return r

# Load meshes and mesh-functions used in the MultiMesh from file
multimesh = MultiMesh()
mfs = []
meshes = []
for i in range(2):
    mesh_i = Mesh()
    with XDMFFile("meshes/multimesh_%d.xdmf" %i) as infile:
        infile.read(mesh_i)
    mvc = MeshValueCollection("size_t", mesh_i, 1)
    with XDMFFile("meshes/mf_%d.xdmf" %i) as infile:
        infile.read(mvc, "name_to_read")
    mfs.append(cpp.mesh.MeshFunctionSizet(mesh_i, mvc))
    meshes.append(mesh_i)
    multimesh.add(mesh_i)

multimesh.build()
multimesh.auto_cover(0,Point(1.25, 0.875))

V = MultiMeshFunctionSpace(multimesh, "CG", 1)
mf_0 = mfs[0]
mf_1 = mfs[1]

x0 = SpatialCoordinate(meshes[0])
x1 = SpatialCoordinate(meshes[1])
T = MultiMeshFunction(V)
T.assign_part(0, project(Constant(0), FunctionSpace(meshes[0], "CG", 1)))
T.assign_part(1, project(x1[0]*x1[1], FunctionSpace(meshes[1], "CG", 1)))

def deformation_vector():
    from femorph import VolumeNormal
    n1 = VolumeNormal(multimesh.part(1))
    bc = DirichletBC(VectorFunctionSpace(multimesh.part(1), "CG",1),
                     Constant((0,0)), mfs[1],2)
    bc.apply(n1.vector())
    S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
    s = MultiMeshFunction(S)
    s.assign_part(1,n1)
    return s

Jsingle = assemble(T.part(1,deepcopy=True)*ds(subdomain_data=mfs[1],
                                              subdomain_id=1))
Jmulti = assemble_multimesh(T("-")*dI)
print(Jsingle, Jmulti)
assert(numpy.isclose(Jsingle,Jmulti))


# Compute gradient
def tan_div(s, n):
    return div(s)-dot(dot(grad(s), n), n)

def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T, n) - dot(grad(s).T, n)
S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
s = TestFunction(S)
n = FacetNormal(multimesh)
dJdOmega = tan_div(s("-"), n("-"))*T("-")*dI

def JT(T):
    return T*dI
J = JT(T("-"))
