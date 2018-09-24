from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import moola
import os
set_log_level(LogLevel.ERROR)
os.system("mkdir -p results")
os.system("mkdir -p figures")

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

# Creating dx/dtheta around center (1.25,0.875)
sx = "-x[1]+0.875"
sy = "x[0]-1.25"
s = Expression((sx,sy),degree=3)

# Source expression
fexp = Expression('x[0]*sin(x[0])*cos(x[1])', degree=4)


def a_s(T, v):
    return dot(grad(T), grad(v))*dX
def l_s(f,v):
    return f*v*dX
def a_N(T,v,n,h, alpha=4.0):
    return - (dot(avg(grad(T)), jump(v, n))*dI 
              + dot(avg(grad(v)), jump(T, n))*dI)+alpha/h*jump(T)*jump(v)*dI
def a_O(T,v, beta=4.0):
    return beta*dot(jump(grad(T)), jump(grad(v)))*dO

def JT(T):
    """
    Returns functional as ufl-expression for given T
    """
    return 0.5*T*T*dX

def eval_J(T):
    """
    Returns functional value with given T
    """
    return 0.5*assemble_multimesh(T**2*dX)

def eval_dJ(T, lmb):
    """
    Computes the shape gradient for the optimization problem given state 
    and adjoint solution
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

    T1 = T.part(1, deepcopy=True)
    lmb1 = lmb.part(1, deepcopy=True)
    dS = Measure("ds", domain=multimesh.part(1), subdomain_data=mf_1)

    normal = FacetNormal(multimesh.part(1))
    d = (0.5*T1*T1-dot(grad(lmb1),normal)*dot(grad(T1),normal))

    # Apply deformation s
    dJs = assemble(inner(normal,s)*d*dS(2))
    return dJs

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
    beta = 4.0

    # Define bilinear form and linear form 
    a = a_s(T,v)+a_N(T,v,n,h)+a_O(T,v)

    f = MultiMeshFunction(V)
    f.interpolate(fexp)
    L = l_s(f,v)

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
    f = MultiMeshFunction(V)
    f.interpolate(fexp)
    n = FacetNormal(multimesh)
    h = 2.0*Circumradius(multimesh)
    h = (h('+') + h('-')) / 2
    v = MultiMeshFunction(V)
    L = JT(T) + a_s(T,v) + a_N(T,v,n,h) + a_O(T,v) - l_s(f,v)
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


def all_angles(fac=1):
    """
    Computes all rotations (360 degree rotation with input fac as
    degree interval)
    """
    T0 = solve_poisson(multimesh)
    lmb0 = solve_adjoint(T0)
    it = 0
    J_old = eval_J(T0)
    i = 0
    degree= 0
    Js = []
    degrees = []
    while i < fac*360:
        multimesh.build()
        T = solve_poisson(multimesh)
        lmb = solve_adjoint(T)
        J = eval_J(T)
        if J < J_old:
            degree = 1./fac*i
            J_old = J
            print(J, 1./fac*i)
        meshes[1].rotate(1./fac, 2, Point(1.25,0.875))
        multimesh.build()
        Js.append(J)
        degrees.append(1./fac*i)
        i+=1
      
    np.savez("results/Global_Dirichlet.npz", deg=degrees, J=Js)
    print(J_old, degree)

def deform_mesh(s, forget=False):
    """
    Rotates the top mesh with angle 's' around (1.25,0.875).
    if forget:
        creates a tmp meshes and computes the functional value
        with the given rotation
    else:
        rotates the multimesh with angle `s`
    """
    if forget:
        mesh_1_new = Mesh(multimesh.part(1))
        multimesh_1 = MultiMesh()
        multimesh_1.add(Mesh(multimesh.part(0)))
        multimesh_1.add(mesh_1_new)
        multimesh_1.build()
    else:
        multimesh_1 = multimesh
        mesh_1_new = meshes[1]

    mesh_1_new.rotate(s,2, Point(1.25,0.875))
    multimesh_1.build()

    if forget:
        T_test = solve_poisson(multimesh_1)
        return eval_J(T_test) 
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    all_angles()
    files = np.load("results/Global_Dirichlet.npz")
    fig = plt.figure()
    plt.plot(files["deg"], files["J"], '-',color=plt.cm.coolwarm(0), linewidth=4)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.xlim(0., 360.0)
    plt.xlabel(r"$\theta$", fontsize=35)
    plt.ylabel(r"$J(T(\theta)))$", fontsize=35,rotation=90)
    plt.grid()
    plt.annotate(r'b)', xy=(1.25, -0.2), color="k",size=20)
    fig.set_size_inches((12, 8.5), forward=False)
    plt.savefig("figures/MultiMeshObstacleAll.png",bbox_inches='tight',format='png', dpi=300)

    os.system("convert figures/MultiMeshObstacleAll.png -trim figures/MultiMeshObstacleAll.png")
    it = 0
    tot_rot = 0
    Js = []
    search = moola.linesearch.ArmijoLineSearch(start_stp=180, stpmin=1e-6)
    while it <100:
        print("Iteration: %d" %it)
        T = solve_poisson(multimesh)
        lmb = solve_adjoint(T)
        J = eval_J(T)
        Js.append(J)
        dJ = eval_dJ(T, lmb) # 1/radians
        print("Functional value: %.3e" % J)
        # #Taylor test
        # def phi_taylor(step):
        #     print("Evaluating functional with degree %.2e" %(step))
        #     J_step = deform_mesh(step, True)
        #     return J_step        
        # def phi_dphi0_taylor():
        #     return float(J), dJ

        # step = 1 # Degrees
        # J0,g = phi_dphi0_taylor()
        # J1 = phi_taylor(step)
        # step = float(step)*2*np.pi/360 # Radians
        # print(g, (J1-J0)/(step))
        # exit(1)
        def phi(step):
            step = - dJ*(180/np.pi)*step
            print("Evaluating functional with degree %.2e" %(step))
            J_step = deform_mesh(step, True)
            return J_step
        def phi_dphi0():
            return float(J), -np.abs(dJ*np.pi/180)

        try:
            step = search.search(phi, None, phi_dphi0())
            step = -dJ*180/np.pi*step # 1/Degree
            print(step)
        except:
            print("Linesearch done")
            break
        
        deform_mesh(step)
        tot_rot += step
        print("Rotation: %.3f" %tot_rot)
        it+=1
    
    print("Total rotation: %.3f" %tot_rot)
    multimesh.build()
    T = solve_poisson(multimesh)
    meshes[1].rotate(-tot_rot)
    multimesh.build()
    meshes[1].rotate(tot_rot)
    multimesh.build()
    #m = np.argmin(files["J"])
    #print(files["deg"][m]-360)
    # plt.plot(Js)

    import matplotlib as mpl
    T_min, T_max = np.min(T.vector().get_local()),\
                   np.max(T.vector().get_local())
    fig,ax = plt.subplots(1,1)
    p=plot(T.part(0,deepcopy=True), norm=mpl.colors.Normalize(vmin=T_min, vmax=T_max),zorder=997)
    p1=plot(T.part(1,deepcopy=True), norm=mpl.colors.Normalize(vmin=T_min, vmax=T_max),zorder=999)

    #Finding center for removal circle
    x = SpatialCoordinate(multimesh.part(1));
    cx = assemble(x[0]*dx(domain=multimesh.part(1)))\
         /assemble(1*dx(domain=multimesh.part(1)))
    cy = assemble(x[1]*dx(domain=multimesh.part(1)))\
         /assemble(1*dx(domain=multimesh.part(1)))
    circ=plt.Circle((cx,cy),0.3,color= mpl.colors.get_named_colors_mapping()['darkgrey'], linewidth=0,alpha=1,zorder=998)
    # ax.add_patch(circ)
    #plt.annotate(r'a)', xy=(1.25, -0.2), color="k",size=20)
    plt.axis([0,2.5, -0.2, 1.75])
    plt.axis("off")
    plt.subplots_adjust(right=0.8)
    cbar_ax = plt.gcf().add_axes([0.81, 0.3, 0.03, 0.45])
    cbar = mpl.colorbar.ColorbarBase(cbar_ax,
                                     norm=mpl.colors.Normalize(vmin=T_min, vmax=T_max))


    plt.savefig("figures/SimpleExampleOptimal.png",bbox_inches='tight',format='png', dpi=300)
    import os
    os.system("convert figures/SimpleExampleOptimal.png -trim figures/SimpleExampleOptimal.png")

