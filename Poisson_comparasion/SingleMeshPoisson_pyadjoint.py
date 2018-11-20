from dolfin import *
from dolfin_adjoint import *
from create_multiple_meshes import inner_marker, outer_marker
import matplotlib.pyplot as plt

# Initialize mesh and facet function
m_name = "meshes/singlemesh.xdmf"
f_name = "meshes/mf.xdmf"

mesh = Mesh()
with XDMFFile(m_name) as infile:
    infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(f_name) as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Move mesh
S = VectorFunctionSpace(mesh, "CG", 1)        
s = Function(S)
ALE.move(mesh, s)

V = FunctionSpace(mesh, "CG", 1)
T = Function(V, name="state")
u = TrialFunction(V)
v = TestFunction(V)
        
# Define bilinear form and linear form
def a_s(T, v):
    return dot(grad(T), grad(v))*dx
def l_s(f, v):
    return f*v*dx
a = a_s(u,v)
x = SpatialCoordinate(mesh)
f = x[0]*sin(x[0])*cos(x[1])
L = l_s(f,v)
dB = Measure("ds", domain=mesh, subdomain_data=mf)
n = FacetNormal(mesh)
alpha = Constant(1)

F = a - l_s(f,v)

A = assemble(lhs(F))
b = assemble(rhs(F))

# Assemble linear system
bcs = [DirichletBC(V, Constant(0), mf, outer_marker),
       DirichletBC(V, Constant(1), mf, inner_marker)]
[bc.apply(A,b) for bc in bcs]
solve(A, T.vector(), b,'lu')

def J_ufl(T):
    """
    Returns functional as ufl-expression for given T
    """
    return 0.5*T*T*dx

J = assemble(J_ufl(T))
Jhat = ReducedFunctional(J, Control(s))
def save_to_file_l2(gradient):
    s.vector()[:] = gradient.get_local()
    File("output/sm_padj_s.pvd") << s
    return s

s1 = Function(S)
s1.vector()[:] = 1
dJ = Jhat.derivative(options={"riesz_representation":save_to_file_l2})

s0 = Function(S)
taylor_test(Jhat, s0, s1,dJdm=s1._ad_dot(dJ))
