from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
import os
from IPython import embed
set_log_level(LogLevel.ERROR)
os.system("mkdir -p results")
os.system("mkdir -p figures")


class PoissonSolver():
    def __init__(self, meshname, facetname, source):
        """
        Initialize Poisson Solver.
        Assumes that
        Arguments:
           meshname (str)    Filename of mesh
           facetname (str)   Filename of facet functions
           source            Source expression or dolfin function
        """
        from create_multiple_meshes import inner_marker, outer_marker
        self.inner_marker = inner_marker
        self.outer_marker = outer_marker
        self.out = File("output/all_track_single.pvd")
        self.f = source
        self.init_mesh(meshname, facetname)
        self.V = FunctionSpace(self.mesh, "CG", 1)
        self.S = VectorFunctionSpace(self.mesh, "CG", 1)        
        self.T = Function(self.V, name="state")
        self.lmb = Function(self.V, anme="adjoint")
        
    def init_mesh(self, meshes_n, facet_funcs): 
        self.mesh = Mesh()
        with XDMFFile(meshes_n) as infile:
            infile.read(self.mesh)
        mvc = MeshValueCollection("size_t", self.mesh, 1)
        with XDMFFile(facet_funcs) as infile:
            infile.read(mvc, "name_to_read")
        self.mf = cpp.mesh.MeshFunctionSizet(self.mesh, mvc)

    def update_mesh(self, perturbation):
        """
        Rotates the top mesh with specified angle (in degrees)
        from initial orientation
        """
        raise NotImplementedError("Single mesh update not implemented")

    # Helper functions for state and adjoint
    def a_s(self, T, v):
        return dot(grad(T), grad(v))*dx
    def l_s(self, f, v):
        return f*v*dx

    def J_ufl(self,T):
        """
        Returns functional as ufl-expression for given T
        """
        return 0.5*T*T*dx

    def eval_J(self, s):
        """
        Evaluates functional with object rotated at given angle
        """
        #self.update_mesh(s)

        # Define trial and test functions and right-hand side
        T = TrialFunction(self.V)
        v = TestFunction(self.V)
        
        # Define bilinear form and linear form 
        a = self.a_s(T,v)
        f = Function(self.V)
        f.interpolate(self.f)
        L = self.l_s(f,v)
        
        # Assemble linear system
        A = assemble(a)
        b = assemble(L)
        bcs = [DirichletBC(self.V, Constant(0), self.mf, self.outer_marker),
               DirichletBC(self.V, Constant(1), self.mf, self.inner_marker)]
        [bc.apply(A,b) for bc in bcs]
        
        # Solving linear system
        solve(A, self.T.vector(), b,'lu')
        self.out << self.T

         # Assemble functional value
        return assemble(self.J_ufl(self.T))


    def eval_dJ(self, angle): # in degrees
        """
        Computes gradient at given angle
        """

        # Update state and mesh
        self.eval_J(angle)

        # Solve adjoint eq
        f = Function(self.V)
        f.interpolate(self.f)
        v = Function(self.V)
        L = self.J_ufl(self.T) + self.a_s(self.T,v) - self.l_s(f,v)
        adjoint = derivative(L, self.T, TestFunction(self.V))
        from ufl import replace
        adjoint = replace(adjoint,  {v: TrialFunction(self.V)})
        a, L = lhs(adjoint), rhs(adjoint)
    
        A = assemble(a)
        b = assemble(L)
        bcs = [DirichletBC(self.V, Constant(0), self.mf, self.outer_marker)
               ,DirichletBC(self.V, Constant(0), self.mf, self.inner_marker)]
        [bc.apply(A,b) for bc in bcs]
        solve(A, self.lmb.vector(), b, 'lu')

        # Compute gradient
        s = TestFunction(self.S)
        d = div(s)*(0.5*self.T*self.T+dot(grad(self.T),grad(self.lmb))
                    -f*self.lmb)*dx - inner(dot(grad(self.T),grad(s)), grad(self.lmb))*dx\
                    - inner(grad(self.T), dot(grad(self.lmb),grad(s)))*dx
        n = FacetNormal(self.mesh)
        d_Hadamard = inner(s,n)*(0.5*self.T*self.T-inner(n, grad(self.lmb))*inner(n, grad(self.T)))*ds

        dJs = assemble(d)
        dJs_Hadamard = assemble(d_Hadamard)
        s = Function(self.S)
        s.vector()[:] = dJs.get_local()
        File("output/sm_s.pvd") << s
        s.vector()[:] = dJs_Hadamard.get_local()
        File("output/sm_h.pvd") << s




    
if __name__ == '__main__':
    m_names = "meshes/singlemesh.xdmf"
    f_names = "meshes/mf.xdmf"
    fexp = Expression('x[0]*sin(x[0])*cos(x[1])', degree=4)
    solver = PoissonSolver(m_names, f_names, fexp)
    solver.eval_dJ(0)

    
