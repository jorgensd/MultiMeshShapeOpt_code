from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import os
from IPython import embed
set_log_level(LogLevel.ERROR)
os.system("mkdir -p results")
os.system("mkdir -p figures")


class PoissonSolver():
    def __init__(self, p, theta, mesh_names, facet_func_names, source):
        """
        Initialize Poisson Solver.
        Assumes that
        Arguments:
           p (dolfin.Point)  Point of rotation
           mesh_names [str, str] List of filenames for back and front mesh
           facet_func_names [str, str] List of filenames for facet functions
           source - Source expression or dolfin function
        """
        self.out = [File("output/all_track.pvd"),File("output/all_track1.pvd")]
        self.f = source
        self.point = p
        self.s = Expression(("-x[1]+%s" % p[1], "x[0]-%s" % p[0]),degree=3)
        self.theta = theta
        self.init_multimesh(mesh_names, facet_func_names, self.theta,
                            self.point)
        self.V = MultiMeshFunctionSpace(self.multimesh, "CG", 1)
        self.T = MultiMeshFunction(self.V, name="state")
        self.lmb = MultiMeshFunction(self.V, anme="adjoint")
        
    def init_multimesh(self, meshes_n, facet_funcs,theta,p): 
        multimesh = MultiMesh()
        mfs = []
        meshes = []
        for i in range(2):
            mesh_i = Mesh()
            with XDMFFile(meshes_n[i]) as infile:
                infile.read(mesh_i)
            mvc = MeshValueCollection("size_t", mesh_i, 1)
            with XDMFFile(facet_funcs[i]) as infile:
                infile.read(mvc, "name_to_read")
            mfs.append(cpp.mesh.MeshFunctionSizet(mesh_i, mvc))
            if i==2:
                mesh_i.rotate(theta, 2, p)
            meshes.append(mesh_i)
            multimesh.add(mesh_i)
        multimesh.build()
        self.mfs = mfs
        self.meshes = meshes
        self.multimesh = multimesh


    def update_mesh(self, angle):
        """
        Rotates the top mesh with specified angle (in degrees)
        from initial orientation
        """
        if not (isinstance(angle, float) or isinstance(angle, int)):
            angle = angle[0]
        
        self.meshes[1].rotate(angle-self.theta, 2, self.point)
        self.multimesh.build()
        self.theta = angle

    # Helper functions for state and adjoint
    def a_s(self, T, v):
        return dot(grad(T), grad(v))*dX
    def l_s(self, f, v):
        return f*v*dX
    def a_N(self, T, v, n, h, alpha=4.0):
        return - (dot(avg(grad(T)), jump(v, n))*dI 
                  + dot(avg(grad(v)), jump(T, n))*dI)\
                  +alpha/h*jump(T)*jump(v)*dI
    def a_O(self, T, v, beta=4.0):
        return beta*dot(jump(grad(T)), jump(grad(v)))*dO

    def J_ufl(self,T):
        """
        Returns functional as ufl-expression for given T
        """
        return 0.5*T*T*dX

    def eval_J(self, angle):
        """
        Evaluates functional with object rotated at given angle
        """
        self.update_mesh(angle)
        mf_0, mf_1 = self.mfs

        # Define trial and test functions and right-hand side
        T = TrialFunction(self.V)
        v = TestFunction(self.V)
    
        # Define facet normal and mesh size
        n = FacetNormal(self.multimesh)
        h = 2.0*Circumradius(self.multimesh)
        h = (h('+') + h('-')) / 2

        # Set parameters
        alpha = 4.0
        beta = 4.0
        
        # Define bilinear form and linear form 
        a = self.a_s(T,v)+self.a_N(T,v,n,h)+self.a_O(T,v)
        
        f = MultiMeshFunction(self.V)
        f.interpolate(self.f)
        L = self.l_s(f,v)
        
        # Deactivate hole in background mesh
        self.multimesh.auto_cover(0, self.point)

        # Assemble linear system
        A = assemble_multimesh(a)
        b = assemble_multimesh(L)
        bcs = [MultiMeshDirichletBC(self.V, Constant(0), mf_0, 1, 0)
               ,MultiMeshDirichletBC(self.V, Constant(1), mf_1, 2, 1)]
        [bc.apply(A,b) for bc in bcs]
        
        # Solving linear system
        self.V.lock_inactive_dofs(A, b)
        solve(A, self.T.vector(), b,'lu')
        self.out[0] << self.T.part(0)
        self.out[1] << self.T.part(1)

         # Assemble functional value
        return assemble_multimesh(self.J_ufl(self.T))


    def eval_dJ(self, angle): # in degrees
        """
        Computes gradient at given angle
        """

        # Update state and mesh
        self.eval_J(angle)

        # Solve adjoint eq
        mf_0, mf_1 = self.mfs
        f = MultiMeshFunction(self.V)
        f.interpolate(self.f)
        n = FacetNormal(self.multimesh)
        h = 2.0*Circumradius(self.multimesh)
        h = (h('+') + h('-')) / 2
        v = MultiMeshFunction(self.V)
        L = self.J_ufl(self.T) + self.a_s(self.T,v)\
            + self.a_N(self.T,v,n,h) + self.a_O(self.T,v) - self.l_s(f,v)
        adjoint = derivative(L, self.T, TestFunction(self.V))
        from ufl import replace
        adjoint = replace(adjoint,  {v: TrialFunction(self.V)})
        a, L = lhs(adjoint), rhs(adjoint)
        self.multimesh.build()
        self.multimesh.auto_cover(0,self.point)
    
        A = assemble_multimesh(a)
        b = assemble_multimesh(L)
        bcs = [MultiMeshDirichletBC(self.V, Constant(0), mf_0, 1, 0)
               ,MultiMeshDirichletBC(self.V, Constant(0), mf_1, 2, 1)]
        [bc.apply(A,b) for bc in bcs]
        self.V.lock_inactive_dofs(A, b)
        solve(A, self.lmb.vector(), b, 'lu')

        # Compute gradient
        T1 = self.T.part(1, deepcopy=True)
        lmb1 = self.lmb.part(1, deepcopy=True)
        dS = Measure("ds", domain=self.multimesh.part(1), subdomain_data=mf_1)
        
        normal = FacetNormal(self.multimesh.part(1))
        d = (0.5*T1*T1-dot(grad(lmb1),normal)*dot(grad(T1),normal))
    
        # Apply deformation s
        dJs = assemble(inner(normal,self.s)*d*dS(2))
        return np.array([180./pi*dJs], dtype=float)



def all_angles():
    import numpy as np
    delta = 5
    N = int(360/delta)
    angles = [delta*i for i in range(N)]
    p = Point(1.25,0.875)
    m_names = ["meshes/multimesh_%d.xdmf" %i for i in range(2)]
    f_names = ["meshes/mf_%d.xdmf" %i for i in range(2)]
    fexp = Expression('x[0]*sin(x[0])*cos(x[1])', degree=4)
    solver = PoissonSolver(p, 0, m_names, f_names, fexp)
    Js = []
    dJds = []
    out0 = File("output/T0_new.pvd")
    out = File("output/T_new.pvd")
    for theta in angles:
        Js.append(solver.eval_J(theta))
        dJds.append(solver.eval_dJ(theta))
        out0 << solver.T.part(0)
        out << solver.T.part(1)
    print(angles[np.argmin(Js)], Js[np.argmin(Js)])
    print(angles[np.argmin(np.abs(dJds))], dJds[np.argmin(np.abs(dJds))])

    np.savez("results/Global_Dirichlet.npz", deg=angles, J=Js)
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
    fig2 = plt.figure()
    plt.plot(angles, dJds, '-',color=plt.cm.coolwarm(0), linewidth=4)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(0., 360.0)
    plt.grid()
    fig.set_size_inches((12, 8.5), forward=False)
    plt.savefig("figures/MultiMeshGrad.png",
                bbox_inches='tight',format='png', dpi=300)

    # os.system("convert figures/MultiMeshObstacleAll.png -trim figures/MultiMeshObstacleAll.png")

    
if __name__ == '__main__':
    all_angles()
    # p = Point(1.25,0.875)
    # m_names = ["meshes/multimesh_%d.xdmf" %i for i in range(2)]
    # f_names = ["meshes/mf_%d.xdmf" %i for i in range(2)]
    # fexp = Expression('x[0]*sin(x[0])*cos(x[1])', degree=4)
    # solver = PoissonSolver(p, 0, m_names, f_names, fexp)


    
