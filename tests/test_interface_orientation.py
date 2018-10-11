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

dI_single = Measure("ds", subdomain_data=mfs[1], subdomain_id=1)
Jsingle = assemble(T.part(1,deepcopy=True)*dI_single)
Jmulti = assemble_multimesh(T("-")*dI)
# Check that interface integral restricted to top mesh is same as boundary
# integral over the same boundary on a single mesh
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
dJ = assemble_multimesh(dJdOmega).get_local()
s_mm = MultiMeshFunction(S)
s_mm.vector()[:] = dJ
s_mm_1 = s_mm.part(1, deepcopy=True)

S_sm = VectorFunctionSpace(multimesh.part(1), "CG", 1)
s_sm = TestFunction(S_sm)
n_sm = FacetNormal(multimesh.part(1))
dJdsingle = tan_div(s_sm, n_sm)*T.part(1,deepcopy=True)*dI_single
dJ_sm = assemble(dJdsingle).get_local()
# Verify tan_div 
assert(numpy.allclose(s_mm_1.vector().get_local(), dJ_sm))

# Verify sign of normal
Jminus = assemble_multimesh(inner(s("-"), n("-"))*dI)
Jminus_1 = MultiMeshFunction(S)
Jminus_1.vector()[:] = Jminus
Jminus_1 = Jminus_1.part(1,deepcopy=True)
Jsingle = assemble(inner(s_sm, n_sm)*dI_single)
assert(numpy.allclose(Jminus_1.vector().get_local(), Jsingle.get_local()))
