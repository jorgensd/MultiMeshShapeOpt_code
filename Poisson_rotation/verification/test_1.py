from dolfin import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace

# Set parameters
alpha = 4.0
beta = 4.0
def JT(T):
    """
    Returns functional as ufl-expression for given T
    """
    return 0.5*T*T*dX


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

V = MultiMeshFunctionSpace(multimesh, "CG", 1)
mf_0 = mfs[0]
mf_1 = mfs[1]

# Define trial and test functions and right-hand side
T = TrialFunction(V)
v = TestFunction(V)

# Define bilinear form and linear form 
a = inner(grad(T), grad(v))*dX
a += inner(T("+"),v("+"))*dI
X = SpatialCoordinate(multimesh)
L = inner(X[0]*X[1], v)*dX

# Deactivate hole in background mesh
def solve_state():
    multimesh.build()
    multimesh.auto_cover(0, Point(1.25, 0.875))
    
    # Assemble linear system
    A = assemble_multimesh(a)
    b = assemble_multimesh(L)
    bcs = [MultiMeshDirichletBC(V, Constant(0), mf_0, 1, 0)
           ,MultiMeshDirichletBC(V, Constant(1), mf_1, 2, 1)]
    [bc.apply(A,b) for bc in bcs]

    # Solving linear system
    T = MultiMeshFunction(V, name="State")
    V.lock_inactive_dofs(A, b)
    solve(A, T.vector(), b,'lu')
    return T

def deformation_vector():
    from femorph import VolumeNormal
    n1 = VolumeNormal(multimesh.part(1))
    S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
    s = MultiMeshFunction(S)
    s.assign_part(1,n1)
    return s

T = solve_state()
J0 = assemble_multimesh(JT(T))

# Solve adjoint eq
adj = TrialFunction(V)
a_adj = inner(grad(v), grad(adj))*dX + v("+")*adj("+")*dI
L_adj = 2*T*v*dX
A_adj = assemble_multimesh(a_adj)
b_adj = assemble_multimesh(L_adj)
bcs_adj = [MultiMeshDirichletBC(V, Constant(0), mf_0, 1, 0),
           MultiMeshDirichletBC(V, Constant(0), mf_1, 2, 1)]
[bc.apply(A_adj,b_adj) for bc in bcs_adj]
adj = MultiMeshFunction(V, name="Adjoint")
V.lock_inactive_dofs(A_adj, b_adj)
solve(A_adj, adj.vector(), b_adj,'lu')

# Compute gradient
def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)
def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)
S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
s = TestFunction(S)
n = FacetNormal(multimesh)
dJdOmega = -inner(dot(nabla_grad(s), nabla_grad(T)), nabla_grad(adj))*dX
dJdOmega+= -inner(nabla_grad(T), dot(nabla_grad(s), nabla_grad(adj)))*dX
dJdOmega+= -dot(grad(X[0]*X[1]), s)* adj*dX
dJdOmega+= -div(s)*(inner(nabla_grad(T), nabla_grad(adj))
                    -X[0]*X[1]*adj+T*T)*dX
dJdOmega+= tan_div(s("+"), n("+"))*T("+")*adj("+")*dI
dJds_ = assemble_multimesh(dJdOmega)


epsilons = [0.01*0.5**i for i in range(5)]
errors = {"0": [],"1": []}
Js = [J0]
for eps in epsilons:
    s_eps = deformation_vector()
    s_eps.vector()[:] *= eps
    dJds = dJds_.inner(s_eps.vector())
    for i in range(2):
        ALE.move(multimesh.part(i), s_eps.part(i))
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
    T_eps = solve_state()
    J_eps = assemble_multimesh(JT(T_eps))
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
rates0 = convergence_rates(errors["0"], epsilons)
rates1 = convergence_rates(errors["1"], epsilons)
print(rates0)
print(rates1)
# Compute gradient and save to file

u = MultiMeshFunction(S)
u.vector()[:] = assemble_multimesh(dJdOmega)
for i in range(2):
    out = XDMFFile("results/dJdO%d.xdmf" %i)
    out.write(u.part(i))
    out.close()

    
