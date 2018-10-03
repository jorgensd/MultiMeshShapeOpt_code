from dolfin import *
from IPython import embed
from pdb import set_trace
import numpy as np
import matplotlib.pyplot as plt

class StokesSolver():

    def __init__(self, points, thetas, mesh_names, facet_func_names):
        """
        Initialize Stokes solver for objects located at "points"
        with orientation "thetas" with meshes from "mesh_names"
        and facet_functions from "facet_func_names".
        Arguments:
            points list(dolfin.Point) - The rotational centers
            thetas list(float)        - The initial orientations
            mesh_names list(str)      - The mesh-filenames
            facet_func_names list(str)- The facet_func filnames
        """
        
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
        for i in range(self.multimesh.num_parts()):
            Vi = FunctionSpace(self.multimesh.part(i), self.V2)
            Pi = FunctionSpace(self.multimesh.part(i), self.S1)
            ui, pi = self.w.part(0, deepcopy=True).split()

            self.u.assign_part(i, interpolate(ui,Vi))
            self.p.assign_part(i, interpolate(pi,Pi))

    def save_state(self):
        for i in range(self.N+1):
            self.outp[i] << u.part(i)
            self.outp[i] << u.part(i)

    
    def update_mesh(self, angles):
        for i in range(1,self.N+1):
            self.meshes[i].rotate(angles[i-1]-self.thetas[i-1],
                                  2, self.points[i-1])
            self.thetas[i-1] = angles[i-1]
        self.multimesh.build()

    
    def eval_J(self, angles):
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
        noslip_value = inflow_value
        obstacle_value = Expression(("0.0", "0.0"), degree=2)
        V = MultiMeshSubSpace(self.VQ, 0)
        Q = MultiMeshSubSpace(self.VQ, 1)
        bc0 = MultiMeshDirichletBC(V, inflow_value,  mf_0, 1, 0)
        bc1 = MultiMeshDirichletBC(V, noslip_value,  mf_0, 3, 0)
        bcs = [bc0, bc1]
        for i in range(1,self.N+1):
            bcs.append(MultiMeshDirichletBC(V, obstacle_value,
                                            mfs[i-1], 2 ,i))
        # Assemble linear system, apply boundary conditions and solve
        A = assemble_multimesh(a)
        b = assemble_multimesh(L)
        [bc.apply(A, b) for bc in bcs]
        self.VQ.lock_inactive_dofs(A, b)
        solve(A, self.w.vector(), b, "mumps")
        self.splitMMF()
        
points = [Point(0.5,0.25), Point(0.75,0.52)]
thetas = [90, 47]
pre = "meshes/"
meshes = [pre+"multimesh_0.xdmf", pre+"multimesh_1.xdmf",pre+"multimesh_1.xdmf"]
mfs = [pre+"mf_0.xdmf", pre+"mf_1.xdmf",pre+"mf_1.xdmf"]
ss = StokesSolver(points, thetas, meshes, mfs)


def plot_mm(multimesh):
    colors = ["r","b","y","c"]
    for i in range(1,multimesh.num_parts()):
        plot(multimesh.part(i),color=colors[i])


ss.eval_J([0,0])
ss.save_state()
File("output/u0.pvd") << ss.u.part(0)
