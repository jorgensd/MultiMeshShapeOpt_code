from dolfin import *
import matplotlib.pyplot as plt
from IPython import embed
from pdb import set_trace
# MultiMesh stability paraameter
beta = 4.0

# Creating dx/dtheta around center (1.25,0.875)
sx = "-x[1]+0.875"
sy = "x[0]-1.25"
s = Expression((sx,sy),degree=3)

def f(X):
    return X[0]*sin(X[0])*cos(X[1])

def a_s(T, v):
    return dot(grad(T), grad(v))*dX
def l_s(f,v):
    return f*v*dX
def a_IP(T,v,n,h, alpha=4.0):
    return - (dot(avg(grad(T)), jump(v, n))*dI 
              + dot(avg(grad(v)), jump(T, n))*dI)+alpha/h*jump(T)*jump(v)*dI
def a_O(T,v, beta=4.0):
    return beta*dot(jump(grad(T)), jump(grad(v)))*dO

def solve_poisson(multimesh):
    """
    Solves the Poisson problem with Dirichlet Boundary conditions on a given
    multimesh
    """
    # Create function space for the temperature
    V = MultiMeshFunctionSpace(multimesh, "CG", 1)
    mfs = []
    for i in range(2):
        mvc = MeshValueCollection("size_t", multimesh.part(i), 1)
        with XDMFFile("meshes/mf_%d.xdmf" %i) as infile:
            infile.read(mvc, "name_to_read")
        mfs.append(cpp.mesh.MeshFunctionSizet(multimesh.part(i), mvc))
    mf_0 = mfs[0]
    mf_1 = mfs[1]

    # Define trial and test functions and right-hand side
    T = TrialFunction(V)
    v = TestFunction(V)

    # Define facet normal and mesh size
    n = FacetNormal(multimesh)
    h = 2.0*Circumradius(multimesh)
    h = (h('+') + h('-')) / 2

    # Set parameters
    alpha = 4.0
    

    # Define bilinear form and linear form 
    a = a_s(T,v)+a_IP(T,v,n,h)+a_O(T,v)
    X = SpatialCoordinate(multimesh)
    f_ = f(X)
    L = l_s(f_,v)

    # Deactivate hole in background mesh
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

def JT(T):
    """
    Returns functional as ufl-expression for given T
    """
    return 0.5*T*T*dX

def solve_adjoint(T):
    """
    Computes the adjoint solution of the Poisson problem given solution T
    """
    multimesh = T.function_space().multimesh()
    mfs = []
    for i in range(2):
        mvc = MeshValueCollection("size_t", multimesh.part(i), 1)
        with XDMFFile("meshes/mf_%d.xdmf" %i) as infile:
            infile.read(mvc, "name_to_read")
        mfs.append(cpp.mesh.MeshFunctionSizet(multimesh.part(i), mvc))
    mf_0 = mfs[0]
    mf_1 = mfs[1]


    V = T.function_space()
    X = SpatialCoordinate(multimesh)
    f_ = f(X)
    n = FacetNormal(multimesh)
    h = 2.0*Circumradius(multimesh)
    h = (h('+') + h('-')) / 2
    v = MultiMeshFunction(V)
    L = JT(T) + a_s(T,v) + a_IP(T,v,n,h) + a_O(T,v) - l_s(f_,v)
    adjoint = derivative(L, T, TestFunction(V))

    from ufl import replace
    adjoint = replace(adjoint,  {v: TrialFunction(V)})
    lmb = MultiMeshFunction(V)
    a, L = lhs(adjoint), rhs(adjoint)
    multimesh.build()
    multimesh.auto_cover(0,Point(1.25, 0.875))
    
    A = assemble_multimesh(a)
    b = assemble_multimesh(L)
    bcs = [MultiMeshDirichletBC(V, Constant(0), mf_0, 1, 0)
           ,MultiMeshDirichletBC(V, Constant(0), mf_1, 2, 1)]
    [bc.apply(A,b) for bc in bcs]
    V.lock_inactive_dofs(A, b)
    solve(A, lmb.vector(), b, 'lu')
    return lmb

def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)

def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)

def compute_gradient(T, lmb, s):
    multimesh = T.function_space().multimesh()
    X = SpatialCoordinate(multimesh)
    h = 2.0*Circumradius(multimesh)
    h = (h('+') + h('-')) / 2
    f_ = f(X)
    df = grad(f_)
    n = FacetNormal(multimesh)
    
    dJOmega = -dot(dot(grad(T), grad(s)), grad(lmb))*dX\
              +div(s)*dot(grad(T), grad(lmb))*dX\
              -dot(grad(T),dot(grad(lmb), grad(s)))*dX\
              +dot(s, df)*lmb*dX + div(s)*f_*lmb*dX
    dJOmega += div(s)*0.5*T*T*dx
    dJdO = -dot(jump(dot(grad(T), grad(s))), jump(grad(lmb)))*dO\
           +div(s)*dot(jump(grad(T)), jump(grad(lmb)))*dO\
           -dot(jump(grad(T)), jump(dot(grad(lmb),grad(s))))*dO
    dJdI = -tan_div(s, n("+"))*dot(n("+"), avg(grad(T)))*jump(lmb)*dI\
           -tan_div(s, n("+"))*dot(n("+"), avg(grad(lmb))*jump(T))*dI\
           -tan_div(s("+"),n("+"))*beta/h*jump(T)*jump(lmb)*dI\
           -dot(dn_mat(s, n("+")), avg(grad(lmb))*jump(T))*dI\
           -dot(dn_mat(s, n("+")), avg(grad(T))*jump(lmb))*dI\
           +dot(n("+"), avg(dot(grad(T), grad(s)))*jump(lmb))*dI\
           +dot(n("+"), avg(dot(grad(lmb), grad(s)))*jump(T))*dI
    
    return dJOmega + dJdO + dJdI 

def deformation_vector(multimesh, step=1):
    from femorph import VolumeNormal
    n2 = VolumeNormal(multimesh.part(1))
    S_sm = VectorFunctionSpace(multimesh.part(1), "CG", 1)
    bc = DirichletBC(S_sm, step*n2, mfs[1], 2)
    us,vs = TrialFunction(S_sm), TestFunction(S_sm)
    a_ = inner(grad(us),grad(vs))*dx + inner(us,vs)*dx
    deformation = Function(S_sm)
    solve(lhs(a_) == rhs(a_), deformation, bcs=bc)
    return deformation

if __name__ == "__main__":
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
    T = solve_poisson(multimesh)
    J0 = assemble_multimesh(JT(T))
    Js = [J0]
    lmb = solve_adjoint(T)
    S = MultiMeshVectorFunctionSpace(multimesh, "CG", 1)
    s_sm = deformation_vector(multimesh, 1)
    s_mm = MultiMeshFunction(S)
    s_mm.assign_part(1,s_sm)
    s = TestFunction(S)
    dJds = assemble_multimesh(compute_gradient(T, lmb, s))
    dJds = dJds.inner(s_mm.vector())

    epsilon = [0.01*0.5**i for i in range(5)]
    for eps in epsilon:
        s_eps = deformation_vector(multimesh, eps)
        ALE.move(multimesh.part(1), s_eps)
        T_eps = solve_poisson(multimesh)
        J_eps = assemble_multimesh(JT(T_eps))
        Js.append(J_eps)
        s_eps.vector()[:] *= -1
        ALE.move(multimesh.part(1), s_eps)

    # Compute gradient and save to file
    # dJ = compute_gradient(T, lmb, s)
    # u = MultiMeshFunction(S)
    # u.vector()[:] = assemble_multimesh(dJ)
    # for i in range(2):
    #     out = XDMFFile("results/dJdO%d.xdmf" %i)
    #     out.write(u.part(i))
    #     out.close()

    
