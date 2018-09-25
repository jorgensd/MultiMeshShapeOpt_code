from dolfin import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace


def f(X):
    return X[0]*sin(X[0])*cos(X[1])

def a_s(T, v):
    return dot(grad(T), grad(v))*dx
def l_s(f,v):
    return f*v*dx(domain=mesh)

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

def solve_adjoint(T):
    """
    Computes the adjoint solution of the Poisson problem given solution T
    """
    mesh = T.function_space().mesh()

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    V = T.function_space()
    X = SpatialCoordinate(mesh)
    f_ = f(X)
    n = FacetNormal(mesh)
    v = Function(V)
    L = JT(T) + a_s(T,v) - l_s(f_,v)
    adjoint = derivative(L, T, TestFunction(V))

    from ufl import replace
    adjoint = replace(adjoint,  {v: TrialFunction(V)})
    lmb = Function(V)
    a, L = lhs(adjoint), rhs(adjoint)
    
    A = assemble(a)
    b = assemble(L)
    bcs = [DirichletBC(V, Constant(0), mf, 1),
           DirichletBC(V, Constant(0), mf, 2)]
    [bc.apply(A,b) for bc in bcs]
    solve(A, lmb.vector(), b, 'lu')
    return lmb

def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)

def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)

def compute_gradient(T, lmb, s):
    mesh = T.function_space().mesh()
    X = SpatialCoordinate(mesh)
    f_ = f(X)
    df = grad(f_)
    n = FacetNormal(mesh)
    
    dJOmega = -dot(dot(grad(T), grad(s)), grad(lmb))*dx\
              +div(s)*dot(grad(T), grad(lmb))*dx\
              -dot(grad(T),dot(grad(lmb), grad(s)))*dx\
              -dot(s, df)*lmb*dx - div(s)*f_*lmb*dx
    dJOmega += div(s)*0.5*T*T*dx
    return dJOmega

def deformation_vector(mesh):
    from femorph import VolumeNormal
    n2 = VolumeNormal(mesh)
    S_sm = VectorFunctionSpace(mesh, "CG", 1)
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

    T = solve_poisson(mesh)
    J0 = assemble(JT(T))
    Js = [J0]
    lmb = solve_adjoint(T)
    S = VectorFunctionSpace(mesh, "CG", 1)
    s_mm = deformation_vector(mesh)

    s = TestFunction(S)
    dJds = assemble(compute_gradient(T, lmb, s))
    grad_out = Function(S)
    grad_out.vector()[:] = dJds
    File("sm_grad.pvd") << grad_out
    dJds = dJds.inner(s_mm.vector())
    epsilons = [0.01*0.5**i for i in range(5)]
    errors = {"0": [],"1": []}
    for eps in epsilons:
        s_eps = deformation_vector(mesh)
        s_eps.vector()[:] *= eps
        ALE.move(mesh, s_eps)
        T_eps = solve_poisson(mesh)
        J_eps = assemble(JT(T_eps))
        Js.append(J_eps)
        errors["0"].append(abs(J_eps-Js[0]))
        errors["1"].append(abs(J_eps-Js[0]-eps*dJds))
        s_eps.vector()[:] *= -1
        ALE.move(mesh, s_eps)
    print(errors["0"])
    print(dJds)
    print(errors["1"])
    rates0 = convergence_rates(errors["0"], epsilons)
    rates1 = convergence_rates(errors["1"], epsilons)
    exit(1)
    print(rates0)
    print(rates1)
    # Compute gradient and save to file
    # dJ = compute_gradient(T, lmb, s)
    # u = MultiMeshFunction(S)
    # u.vector()[:] = assemble_multimesh(dJ)
    # for i in range(2):
    #     out = XDMFFile("results/dJdO%d.xdmf" %i)
    #     out.write(u.part(i))
    #     out.close()

    
