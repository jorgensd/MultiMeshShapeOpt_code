from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import os
from IPython import embed
set_log_level(LogLevel.ERROR)
os.system("mkdir -p results")
os.system("mkdir -p figures")
from create_multiple_meshes import inner_marker,outer_marker,c_x,c_y

class PoissonSolver():
    def __init__(self, mesh_names, facet_func_names, source):
        """
        Initialize Poisson Solver.
        Assumes that
        Arguments:
           mesh_names [str, str] List of filenames for back and front mesh
           facet_func_names [str, str] List of filenames for facet functions
           source - Source expression or dolfin function
        """
        self.out = [File("output/all_track.pvd"),File("output/all_track1.pvd")]
        self.f = source
        self.init_multimesh(mesh_names, facet_func_names)
        self.V = MultiMeshFunctionSpace(self.multimesh, "CG", 2)
        self.T = MultiMeshFunction(self.V, name="state")
        self.lmb = MultiMeshFunction(self.V, name="adjoint")
        self.S = MultiMeshVectorFunctionSpace(self.multimesh, "CG", 1)
        self.s = MultiMeshFunction(self.S)
        self.point = Point(c_x,c_y)
        
    def init_multimesh(self, meshes_n, facet_funcs): 
        multimesh = MultiMesh()
        mfs = []
        meshes = []
        for i in range(len(meshes_n)):
            mesh_i = Mesh()
            with XDMFFile(meshes_n[i]) as infile:
                infile.read(mesh_i)
            mvc = MeshValueCollection("size_t", mesh_i, 1)
            with XDMFFile(facet_funcs[i]) as infile:
                infile.read(mvc, "name_to_read")
            mfs.append(cpp.mesh.MeshFunctionSizet(mesh_i, mvc))
            meshes.append(mesh_i)
            multimesh.add(mesh_i)
        multimesh.build()
        self.mfs = mfs
        self.meshes = meshes
        self.multimesh = multimesh


    def update_mesh(self, s):
        """
        Moves top mesh in direction s
        """
        ALE.move(self.multimesh.part(1), s)
        self.multimesh.build()

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

    def eval_J(self, s):
        """
        Evaluates functional with object rotated at given angle
        """
        self.update_mesh(s)
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


    def eval_dJ(self, s): # in degrees
        """
        Computes gradient at given angle
        """

        # Update state and mesh
        self.eval_J(s)

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

        # Add outer boundary_terms
        dJ = 0
        for i in range(1,self.multimesh.num_parts()):
            Ti = self.T.part(i, deepcopy=True)
            lmbi = self.lmb.part(i, deepcopy=True)
            ni = FacetNormal(self.multimesh.part(i))
            si = TestFunction(self.s.part(i,deepcopy=True).function_space())
            dsi = Measure("ds", domain=self.multimesh.part(i),
                         subdomain_data=mf_1)
            dJ += assemble(inner(ni,si)*(Ti**2-dot(grad(lmbi), ni)*dot(grad(Ti),ni))*dsi(inner_marker))
        s1 = self.s.part(1,deepcopy=True)
        s1.vector()[:] = dJ.get_local()

        s = TestFunction(self.S)
        n = FacetNormal(self.multimesh)
        # n("+") is n on background mesh, n("-") on top mesh
        def tan_div(s, n):
            return div(s)-dot(dot(grad(s),n),n)
        dJ_full = assemble_multimesh(
            inner(n("-"),s("-"))*(
                avg(self.lmb)*jump(-div(grad(self.T)))
                +self.lmb("-")*dot(dot(jump(grad(grad(self.T))),n("-")),n("-"))
                -tan_div(jump(grad(self.T))*self.lmb("-"), n("-")))
                *dI)
        print(np.max(s1.vector().get_local()))
        print(np.max(dJ_full.get_local()))
        s_tmp = MultiMeshFunction(self.S)
        s_tmp.vector()[:] = dJ_full.get_local()
        s1 += s_tmp.part(1,deepcopy=True)
        return s1



    
if __name__ == '__main__':
    m_names = ["meshes/multimesh_%d.xdmf" %i for i in range(2)]
    f_names = ["meshes/mf_%d.xdmf" %i for i in range(2)]
    fexp = Expression('x[0]*sin(x[0])*cos(x[1])', degree=4)
    solver = PoissonSolver(m_names, f_names, fexp)
    s_org = solver.s.part(1,deepcopy=True)
    dJ = solver.eval_dJ(s_org)
    from IPython import embed; embed()

    
