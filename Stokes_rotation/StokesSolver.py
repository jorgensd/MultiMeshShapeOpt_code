from dolfin import *
from IPython import embed
from pdb import set_trace
import numpy as np
import matplotlib.pyplot as plt

class StokesSolver():
    set_log_level(LogLevel.ERROR)
    def __init__(self, points, thetas, mesh_names, facet_func_names, inlets):
        """
        Initialize Stokes solver for objects located at "points"
        with orientation "thetas" with meshes from "mesh_names"
        and facet_functions from "facet_func_names".
        Arguments:
            points list(dolfin.Point) - The rotational centers
            thetas list(float)        - The initial orientations
            mesh_names list(str)      - The mesh-filenames
            facet_func_names list(str)- The facet_func filnames
            inlets dict               - Dictonary containing positions of inlets
                                        as well as amplitude
        """
        self.inlets = inlets
        self.outlet_marker = 3
        self.obstacle_marker = 4
        self.wall_marker = 5
        self.N = len(points)
        self.outu = [File("output/u_%d.pvd" %i) for i in range(self.N+1)]
        self.outp = [File("output/p_%d.pvd" %i) for i in range(self.N+1)]
        self.points = points
        self.s = [Expression(("-x[1]+%s" % points[i][1],
                              "x[0]-%s" % points[i][0]), degree=3)
                  for i in range(self.N)]
        self.thetas = thetas
        self.init_multimesh(mesh_names, facet_func_names, self.thetas,
                            self.points)
        self.V2 = VectorElement("CG", triangle, 2)
        self.S1 = FiniteElement("CG", triangle, 1)
        self.VQ = MultiMeshFunctionSpace(self.multimesh, self.V2*self.S1)
        V = MultiMeshFunctionSpace(self.multimesh, self.V2)
        Q = MultiMeshFunctionSpace(self.multimesh, self.S1)
        self.w = MultiMeshFunction(self.VQ, name="State")
        self.u = MultiMeshFunction(V, name="u")
        self.p = MultiMeshFunction(Q, name="p")

    def init_multimesh(self, meshes_n, facet_funcs,theta,p): 
        multimesh = MultiMesh()
        mfs = []
        meshes = []
        for i in range(self.N+1):
            mesh_i = Mesh()
            with XDMFFile(meshes_n[i]) as infile:
                infile.read(mesh_i)
            mvc = MeshValueCollection("size_t", mesh_i, 1)
            with XDMFFile(facet_funcs[i]) as infile:
                infile.read(mvc, "name_to_read")
            mfs.append(cpp.mesh.MeshFunctionSizet(mesh_i, mvc))
            if i>0:
                mesh_i.translate(Point(p[i-1][0]-0.5, p[i-1][1]-0.5))
                mesh_i.rotate(theta[i-1], 2, p[i-1])
            meshes.append(mesh_i)
            multimesh.add(mesh_i)
        multimesh.build()
        self.mfs = mfs
        self.meshes = meshes
        self.multimesh = multimesh

    """
    Several helper functions for the linear and bilinear weak formulation
    of the Stokes equation.
    (see https://arxiv.org/pdf/1206.1933.pdf p. 5 eq 3.1-3.2)
    """
    def b_h(self, v, q, n):
        return -div(v)*q*dX + jump(v, n)*avg(q)*dI
    
    def l_h(self, v, q, f):
        return inner(f, v)*dX

    def s_O(self, v, w):
        return inner(jump(grad(v)), jump(grad(w)))*dO
    
    def s_C(self, v, q, w, r, h):
        return h*h*inner(-div(grad(v)) + grad(q), -div(grad(w))
                         - grad(r))*dC + \
                         h("+")*h("+")*inner(-div(grad(v("+")))
                                             + grad(q("+")),
                                             -div(grad(w("+")))
                                             - grad(r("+")))*dO
    def l_C(self,v, q, f, h):
        return h*h*inner(f, -div(grad(v)) - grad(q))*dC\
            +  h("+")*h("+")*inner(f("+"), -div(grad(v("+")))
                                   - grad(q("+")))*dO
    
    def tensor_jump(self, v, n):
        return outer(v('+'), n('+')) + outer(v('-'), n('-'))
    
    def a_h(self,v, w, n, h,  alpha=6.0): 
        return inner(grad(v), grad(w))*dX \
            - inner(avg(grad(v)), self.tensor_jump(w, n))*dI \
            - inner(avg(grad(w)), self.tensor_jump(v, n))*dI \
            + alpha/avg(h) * inner(jump(v), jump(w))*dI

    def splitMMF(self):
        """
        Split a mixed multimeshfunction into separate multimeshfunctions
        """
        for i in range(self.multimesh.num_parts()):
            Vi = FunctionSpace(self.multimesh.part(i), self.V2)
            Pi = FunctionSpace(self.multimesh.part(i), self.S1)
            ui, pi = self.w.part(i, deepcopy=True).split()

            self.u.assign_part(i, interpolate(ui,Vi))
            self.p.assign_part(i, interpolate(pi,Pi))

    def save_state(self):
        """
        Save current velocity and pressure to file
        """
        for i in range(self.multimesh.num_parts()):
            self.outu[i] << self.u.part(i)
            self.outp[i] << self.p.part(i)
    
    def update_mesh(self, angles):
        """
        Rotate obstacles to angle specified in arrayx
        """
        for i in range(1,self.N+1):
            self.meshes[i].rotate(angles[i-1]-self.thetas[i-1],
                                  2, self.points[i-1])
            self.thetas[i-1] = angles[i-1]
        self.multimesh.build()

    def ufl_J(self, u):
        return inner(grad(u),grad(u))*dX
    
    def eval_J(self, angles):
        self.multimesh.build()
        self.update_mesh(angles)
        mf_0 = self.mfs[0]
        mfs = self.mfs[1:]
        (u, p) = TrialFunctions(self.VQ)
        (v, q) = TestFunctions(self.VQ)
        f = Constant((0.0, 0.0))
        # Define facet normal and mesh size
        n = FacetNormal(self.multimesh)
        h = 2.0*Circumradius(self.multimesh)
        
        # Define bilinear and linear form
        a = self.a_h(u, v, n, h) + self.b_h(v, p, n) + self.b_h(u, q, n)\
            + self.s_O(u, v) + self.s_C(u, p, v, q, h)
        L  = self.l_h(v, q, f) + self.l_C(v, q, f, h)
        
        # Set inactive dofs
        for i in range(self.N):
            self.multimesh.auto_cover(0, self.points[i])

        # Create boundary conditions
        inflow_value = Expression(("1.0", "0.0"),degree=1)
        outflow_value = Constant(0)
        noslip_value =  Expression(("0.0", "0.0"), degree=2)
        obstacle_value = Expression(("0.0", "0.0"), degree=2)
        V = MultiMeshSubSpace(self.VQ, 0)
        Q = MultiMeshSubSpace(self.VQ, 1)
        bcs = []
        for inlet in self.inlets:
            bcs.append(MultiMeshDirichletBC(V, inlet[0],  mf_0, inlet[1], 0))
        bcs.append(MultiMeshDirichletBC(V, noslip_value,  mf_0,
                                        self.wall_marker, 0))

        for i in range(1,self.N+1):
            bcs.append(MultiMeshDirichletBC(V, obstacle_value,
                                            mfs[i-1], self.obstacle_marker ,i))
        # Assemble linear system, apply boundary conditions and solve
        A = assemble_multimesh(a)
        b = assemble_multimesh(L)
        [bc.apply(A, b) for bc in bcs]
        self.VQ.lock_inactive_dofs(A, b)
        solve(A, self.w.vector(), b, "mumps")
        self.splitMMF()

        return assemble_multimesh(self.ufl_J(self.u))

    def eval_dJ(self,angles):

        self.eval_J(angles)

        dJ = np.zeros(self.N)
        for i in range(1,self.N+1):
            u_i = self.u.part(i, deepcopy=True)
            stokes_i = inner(grad(u_i), grad(u_i))
            mf_i = self.mfs[i]
            normal_i = FacetNormal(self.multimesh.part(i))
            dS_i = Measure("ds", domain=self.multimesh.part(i),
                           subdomain_data=mf_i)
            dJ[i-1] = -assemble(stokes_i*inner(normal_i, self.s[i-1])
                               *dS_i(self.obstacle_marker))
        return 180./pi*dJ


def all_angles():
    # solve single rotation problem for all angles and compute gradient
    import numpy as np
    import os
    os.system("mkdir -p figures")
    angles = np.linspace(0,180,25)
    points = [Point(0.5,0.5)]
    thetas = [0]
    inlet_str= "-A*(x[1]-x_l)*(x[1]-x_u)"
    inlet_data = [[Expression((inlet_str, "0"), x_l=0.1, x_u=0.4,
                                     A=250, degree=5), 1],
                         [Expression((inlet_str, "0"), x_l=0.7, x_u=0.85,
                                     A=0, degree=5), 2]]
    pre = "meshes/"
    meshes = [pre+"multimesh_0.xdmf"] +  [pre+"multimesh_1.xdmf"]*len(points)
    mfs = [pre+"mf_0.xdmf"] + [pre+"mf_1.xdmf"]*len(points)
    solver = StokesSolver(points, thetas, meshes, mfs, inlet_data)

    Js = []
    dJds = []
    out0 = File("output/Stokes0_new.pvd")
    out = File("output/Stokes_new.pvd")
    minmax_theta = []
    minmax_grad = []
    i = 0
    for theta in angles:
        Js.append(solver.eval_J([theta]))
        dJds.append(solver.eval_dJ([theta]))
        if i >0 and dJds[-1]/dJds[-2] < 0:
            print("Minmax found")
            minmax_theta.append(angles[i-1])
            minmax_theta.append(theta)
            minmax_grad.append(dJds[-2])
            minmax_grad.append(dJds[-1])
            print(theta, Js[-1], dJds[-1],dJds[-2])
        out0 << solver.u.part(0)
        out << solver.u.part(1)
        i+=1
    print(angles[np.argmin(Js)], Js[np.argmin(Js)])
    print(angles[np.argmin(np.abs(dJds))], dJds[np.argmin(np.abs(dJds))])

    fig = plt.figure()
    plt.plot(angles, Js, '-',color=plt.cm.coolwarm(0), linewidth=4)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.xlim(min(angles), max(angles))
    plt.xlabel(r"$\theta$", fontsize=35)
    plt.ylabel(r"$J(T(\theta)))$", fontsize=35,rotation=90)
    plt.grid()
    fig.set_size_inches((12, 8.5), forward=False)
    plt.savefig("figures/MultiMeshObstacleAll.png",bbox_inches='tight',format='png', dpi=300)
    fig2 = plt.figure()
    plt.plot(angles, dJds, '-',color=plt.cm.coolwarm(0), linewidth=4)
    plt.plot(minmax_theta, minmax_grad,"ro")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(min(angles), max(angles))
    plt.grid()
    fig.set_size_inches((12, 8.5), forward=False)
    plt.savefig("figures/MultiMeshGrad.png",
                bbox_inches='tight',format='png', dpi=300)


if __name__ == "__main__":
    all_angles()
    exit(1)
    points = [Point(0.5,0.5)]#[Point(0.5,0.25), Point(0.75,0.52), Point(0.3,0.8)]
    thetas = [0]#[90, 47, 32]
    inlet_str= "-A*(x[1]-x_l)*(x[1]-x_u)"
    inlet_data = [[Expression((inlet_str, "0"), x_l=0.1, x_u=0.4,
                                     A=250, degree=5), 1],
                         [Expression((inlet_str, "0"), x_l=0.7, x_u=0.85,
                                     A=0, degree=5), 2]]
    # points = [Point(0.5,0.5)]
    # thetas = [63]

    pre = "meshes/"
    meshes = [pre+"multimesh_0.xdmf"] +  [pre+"multimesh_1.xdmf"]*len(points)
    mfs = [pre+"mf_0.xdmf"] + [pre+"mf_1.xdmf"]*len(points)
    ss = StokesSolver(points, thetas, meshes, mfs, inlet_data)

    
    def plot_mm(multimesh):
        colors = ["r","b","y","c"]
        for i in range(1,multimesh.num_parts()):
            plot(multimesh.part(i),color=colors[i])
        

    ss.eval_dJ(thetas)
    embed()
    ss.save_state()
