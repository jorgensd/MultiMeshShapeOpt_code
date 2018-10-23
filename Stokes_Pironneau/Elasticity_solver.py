from dolfin import *
from IPython import embed
from matplotlib.pyplot import show

def epsilon(u):
    """
    Helper function for variational form
    """
    return sym(grad(u))


class ElasticitySolver():
    """
    Linear elasticity problem for a pure Neumann problem
    See paper by Kutcha for information about this: 
    https://arxiv.org/pdf/1609.09425.pdf
    """
    def __init__(self, mesh, facet_function):
        parameters["linear_algebra_backend"] = "PETSc"
        self.mesh = mesh
        self.mf = facet_function
        self.f = Constant((0,0))
        self.h = Constant((0,0))
        self.V= VectorFunctionSpace(mesh, "CG", 1)
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)
        self.u_ = Function(self.V)
        self.mu = 500
        self.lmb = 0
        self._sigma(self.mu, self.lmb)
        self.build_nullspace()
    
    def build_nullspace(self):
        """Function to build null space for 2D elasticity"""
        x = self.u_.vector()
        
        # Create list of vectors for null space
        nullspace_basis = [x.copy() for i in range(3)]

        # Build translational null space basis
        self.V.sub(0).dofmap().set(nullspace_basis[0], 1.0);
        self.V.sub(1).dofmap().set(nullspace_basis[1], 1.0);
        
        # Build rotational null space basis
        self.V.sub(0).set_x(nullspace_basis[2], -1.0, 1);
        self.V.sub(1).set_x(nullspace_basis[2],  1.0, 0);
        
        for x in nullspace_basis:
            x.apply("insert")

        # Create vector space basis and orthogonalize
        basis = VectorSpaceBasis(nullspace_basis)
        basis.orthonormalize()
        self.null_space = basis

    def set_volume_forces(self, f):
        """
        Set volume forces
        """
        self.f = f
    
    def set_boundary_stress(self, h, marker=None):
        """
        Set boundary stress function with given marker
        """
        self.h = h
        if marker is None:
            self.dstress = ds
        else:
            self.dstress =  Measure("ds", domain=self.mesh,
                                    subdomain_data=self.mf,
                                    subdomain_id=marker)

            
    def _sigma(self, mu, lmb):
        """
        Helper function for  variational form, inititalized with 
        """
        self.sigma = 2*mu*epsilon(self.u) + lmb*tr(epsilon(self.u))*Identity(2)

    def solve(self, f, h, marker):
        self.set_volume_forces(f)
        self.set_boundary_stress(h, marker)
        a = inner(self.sigma, grad(self.v))*dx
        L = inner(self.f,self.v)*dx + inner(self.h,self.v)*self.dstress

        # Assemble system
        A = assemble(a)
        b = Vector(MPI.comm_world, self.V.dim())
        assemble(L, tensor=b)

        # Associate null space with A
        as_backend_type(A).set_nullspace(self.null_space)
        as_backend_type(A).set_near_nullspace(self.null_space)

        # Orthogonalize right-hand side to make sure that input is in the
        # range of A aka the orthogonal complement of the null space, cf.
        # linear algebra 101.
        self.null_space.orthogonalize(b)

        # Define solver, set operator and set nullspace
        solver = PETScKrylovSolver("cg", "hypre_amg")
        solver.parameters["monitor_convergence"] = 1==0
        solver.set_operator(A)
        solver.solve(self.u_.vector(), b)
        # plot(self.u_)
        # show()
        
if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("meshes/multimesh_1.xdmf") as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf_1.xdmf") as infile:
        infile.read(mvc, "name_to_read")
        mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    
    e_solve = ElasticitySolver(mesh, mf)
    h = Expression(("x[0]","x[1]"), degree=1)
    f = Constant(("0","0"))
    e_solve.solve(f,h,2)

