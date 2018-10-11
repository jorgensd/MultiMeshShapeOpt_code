from dolfin import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace
import numpy
# Verification of functional only consisting of an interface term
# that does not couple the PDEs and lives on the top mesh

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
T.assign_part(0, project(sin(x0[1]*x0[0]), FunctionSpace(meshes[0], "CG", 1)))
T.assign_part(1, project(cos(x1[1])*10*x1[1], FunctionSpace(meshes[1], "CG", 1)))
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


# Compute gradient
def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)
def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)
S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
s = TestFunction(S)
n = FacetNormal(multimesh)
def JT(T):
    return T("+")*T("-")*dI
J = JT(T)
dJdOmega = (tan_div(s("-"), n("-"))*T("+")*T("-") # Standard part of derivative
            + T("-")*dot(s("-"), grad(T("+"))))*dI # Material derivative of bottom values


dJds_ = assemble_multimesh(dJdOmega)

epsilons = [0.01*0.5**i for i in range(5)]
errors = {"0": [],"1": []}
Js = [assemble_multimesh(J)]
for eps in epsilons:
    s_eps = deformation_vector()
    s_eps.vector()[:] *= eps
    dJds = dJds_.inner(s_eps.vector())
    for i in range(2):
        ALE.move(multimesh.part(i), s_eps.part(i))
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
    J_eps = assemble_multimesh(JT(T))
    Js.append(J_eps)
    errors["0"].append(abs(J_eps-Js[0]))
    errors["1"].append(abs(J_eps-Js[0]-dJds))
    s_eps.vector()[:] *= -1
    for i in range(2):
        ALE.move(multimesh.part(i), s_eps.part(i))
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
print(errors["0"])
print(errors["1"])
rates0 = convergence_rates(errors["0"], epsilons)
rates1 = convergence_rates(errors["1"], epsilons)
assert(min(rates0)>0.95)
assert(min(rates1)>1.85)
print(rates0)
print(rates1)
