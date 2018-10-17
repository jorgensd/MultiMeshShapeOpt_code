from matplotlib.pyplot import show
from pdb import set_trace
from IPython import embed
from numpy import log
from femorph import VolumeNormal
from dolfin import (assemble_multimesh,
                    ALE,
                    Constant,
                    cpp,
                    DirichletBC,
                    FunctionSpace,
                    Mesh, MeshValueCollection,
                    MultiMesh, MultiMeshFunction, MultiMeshFunctionSpace,
                    MultiMeshVectorFunctionSpace,
                    Point, project, SpatialCoordinate,
                    VectorFunctionSpace,
                    sin, cos,
                    div, dot, grad, outer, inner, nabla_grad,
                    dX, dI, dO,
                    TestFunction, TrialFunction,
                    solve,
                    XDMFFile)

# Helper functions for shape calculus
def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)
def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)


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

def deformation_vector():
    # Defines the deformation of the top mesh
    n1 = VolumeNormal(multimesh.part(1))
    bc = DirichletBC(VectorFunctionSpace(multimesh.part(1), "CG",1),
                     Constant((0,0)), mfs[1],2)
    bc.apply(n1.vector())
    S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
    s = MultiMeshFunction(S)
    s.assign_part(1,n1)
    return s

def project_to_background(s_top):
    # Projects a function living on the boundary intersecting with the
    # bottom mesh to the bottom mesh
    S = s_top.function_space()
    u, v = TrialFunction(S), TestFunction(S)
    a = inner(u("+"), v("+"))*dI
    A = assemble_multimesh(a)
    A.ident_zeros()
    l = inner(s_top("-"), v("+"))*dI
    L = assemble_multimesh(l)
    S.lock_inactive_dofs(A,L)
    s_back = MultiMeshFunction(S)
    solve(A, s_back.vector(), L)
    return s_back

s_top = deformation_vector()
s_bottom = project_to_background(s_top)

# Create Spatial Coordinates for each mesh
x0 = SpatialCoordinate(meshes[0])
x1 = SpatialCoordinate(meshes[1])

degree=2
V = MultiMeshFunctionSpace(multimesh, "CG", degree)
# Assign a function to each mesh, such that T and lmb are discontinuous at the
# interface Gamma
T = MultiMeshFunction(V)
T.assign_part(0, project(sin(x0[1]), FunctionSpace(meshes[0], "CG", degree)))
T.assign_part(1, project(cos(x1[0])*x1[1], FunctionSpace(meshes[1], "CG", degree)))
lmb = MultiMeshFunction(V)
lmb.assign_part(0, project(cos(x0[1])*x0[0]*x0[0],
                           FunctionSpace(meshes[0], "CG", degree)))
lmb.assign_part(1, project(x1[0]*sin(x1[1]),
                           FunctionSpace(meshes[1], "CG", degree)))


# Create bilinear form and corresponding gradients

a1 = inner(grad(T), grad(lmb))*dX
# Classic shape derivative term top mesh
da1_top =  div(s_top)*inner(grad(T), grad(lmb))*dX
# Term stemming from grad(T)
da1_top -= inner(dot(grad(s_top), grad(T)), grad(lmb))*dX
# Term stemming from grad(lmb)
da1_top -= inner(grad(T), dot(grad(s_top), grad(lmb)))*dX
# Classic shape derivative term bottom mesh
da1_bottom =  div(s_bottom)*inner(grad(T), grad(lmb))*dX
# Term stemming from grad(T)
da1_bottom += inner(grad(dot(s_bottom, grad(T))), grad(lmb))*dX
# Term stemming from grad(lmb)
da1_bottom += inner(grad(T), grad(dot(s_bottom, grad(lmb))))*dX
# Material derivative of background T
da1_bottom -= inner(dot(nabla_grad(s_bottom), grad(T)), grad(lmb))*dX
# Material derivative of background lmb
da1_bottom -= inner(grad(T), dot(nabla_grad(s_bottom), grad(lmb)))*dX

dJds = assemble_multimesh(da1_top + da1_bottom)
J = a1


# Do a taylor test for deformation of the top mesh
Js = [assemble_multimesh(J)]
epsilons = [0.001*0.5**i for i in range(5)]
errors = {"0": [], "1": []}
for eps in epsilons:
    s_eps = deformation_vector()
    s_eps.vector()[:] *= eps
    for i in range(2):
        ALE.move(multimesh.part(i), s_eps.part(i))
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
    J_eps = assemble_multimesh(a1)
    Js.append(J_eps)
    errors["0"].append(abs(J_eps-Js[0]))
    errors["1"].append(abs(J_eps-Js[0]-eps*dJds))
    s_eps.vector()[:] *= -1
    for i in range(2):
        ALE.move(multimesh.part(i), s_eps.part(i))
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
print(errors["0"])
print(errors["1"])

def convergence_rates(E_values, eps_values):
    r = []
    for i in range(1, len(eps_values)):
        r.append(log(E_values[i]/E_values[i-1])/log(eps_values[i]/
                                                    eps_values[i-1]))

    print("Computed convergence rates: {}".format(r))
    return r


rates0 = convergence_rates(errors["0"], epsilons)
rates1 = convergence_rates(errors["1"], epsilons)
assert(sum(rates0)/len(rates0)>0.95)
assert(sum(rates1)/len(rates1)>1.95)
print(rates0)
print(rates1)
print(sum(rates1)/len(rates1))
