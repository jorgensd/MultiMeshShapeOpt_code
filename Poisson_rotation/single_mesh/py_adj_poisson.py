from dolfin import *
from dolfin_adjoint import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace


def f(X):
    return X[0]*sin(X[0])*cos(X[1])

def a_s(T, v):
    return dot(grad(T), grad(v))*dx
def l_s(f,v):
    return f*v*dx

def solve_poisson(mesh):
    """
    Solves the Poisson problem with Dirichlet Boundary conditions on a given
    multimesh
    """
    # Create function space for the temperature
    V = FunctionSpace(mesh, "CG", 1)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    # Define trial and test functions and right-hand side
    T = TrialFunction(V)
    v = TestFunction(V)

    # Define facet normal and mesh size
    n = FacetNormal(mesh)
    

    # Define bilinear form and linear form 
    a = a_s(T,v)
    X = SpatialCoordinate(mesh)
    f_ = f(X)
    L = l_s(f_,v)
    
    # Assemble linear system
    A = assemble(a)
    b = assemble(L)
    bcs = [DirichletBC(V, Constant(0), mf, 1),
           DirichletBC(V, Constant(1), mf, 2)]
    [bc.apply(A,b) for bc in bcs]

    # Solving linear system
    T = Function(V, name="State")
    solve(A, T.vector(), b,'lu')
    return T

def JT(T):
    """
    Returns functional as ufl-expression for given T
    """
    return 0.5*T*T*dx

def deformation_vector(mesh, S_sm):
    from femorph import VolumeNormal
    n2 = VolumeNormal(mesh)
    # S_sm = VectorFunctionSpace(mesh, "CG", 1)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    bc = DirichletBC(S_sm, n2, mf, 2)
    bc2 = DirichletBC(S_sm, Constant((0,0)), mf, 1)

    us,vs = TrialFunction(S_sm), TestFunction(S_sm)
    a_ = inner(grad(us),grad(vs))*dx + inner(us,vs)*dx
    deformation = Function(S_sm)
    solve(lhs(a_) == rhs(a_), deformation, bcs=[bc,bc2])
    return deformation

def convergence_rates(E_values, eps_values):
    from numpy import log
    r = []
    for i in range(1, len(eps_values)):
        r.append(log(E_values[i]/E_values[i-1])/log(eps_values[i]/
                                                    eps_values[i-1]))

    print("Computed convergence rates: {}".format(r))
    return r

if __name__ == "__main__":
    # Load meshes and mesh-functions used in the MultiMesh from file

    mesh = Mesh()
    with XDMFFile("meshes/singlemesh.xdmf") as infile:
        infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    S = VectorFunctionSpace(mesh, "CG", 1)
    s = Function(S)
    ALE.move(mesh, s)
    T = solve_poisson(mesh)
    J = assemble(JT(T))
    rf = ReducedFunctional(J, Control(s))
    tape.optimize()
    tape.visualise("tape_pa.dot", dot=True)
    s_mm = deformation_vector(mesh, S)
    result = taylor_to_dict(rf, s, s_mm)
    
    print(result["FD"]["Residual"])
    print(result["dJdm"]["Residual"])

    
