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
    def __init__(self, mesh, facet_function, free_marker, deform_marker, constant_mu=True):
        parameters["linear_algebra_backend"] = "PETSc"
        self.mesh = mesh
        self.mf = facet_function
        self.free_marker = free_marker
        self.deform_marker = deform_marker

        self.f = Constant((0,0))
        self.h = Constant((0,0))
        self.V= VectorFunctionSpace(mesh, "CG", 1)
        self.u = TrialFunction(self.V)
        self.v = TestFunction(self.V)
        self.u_ = Function(self.V)
        self.lmb = 0
        self.compute_mu(constant_mu)
        
        self._sigma()
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

    def compute_mu(self, constant):
        """
        Compute mu as according to arxiv paper
        https://arxiv.org/pdf/1509.08601.pdf
        """
        mu_min=Constant(400)
        mu_max=Constant(500)
        if constant:
                self.mu = mu_max
        else:
            V = FunctionSpace(self.mesh, "CG",1)
            u, v = TrialFunction(V), TestFunction(V)
            a = inner(grad(u),grad(v))*dx
            l = Constant(0)*v*dx
            bcs = [DirichletBC(V, mu_min, self.mf, self.free_marker),
                   DirichletBC(V, mu_max, self.mf, self.deform_marker)]
            mu = Function(V)
            solve(a==l, mu, bcs=bcs)
            File("output/mu.pvd") << mu
            self.mu = mu
        
            
    def _sigma(self):
        """
        Helper function for  variational form, inititalized with 
        """
        self.sigma = 2*self.mu*epsilon(self.u) \
                     + self.lmb*tr(epsilon(self.u))*Identity(2)

    def solve(self, f, h):
        self.set_volume_forces(f)
        self.set_boundary_stress(h, self.deform_marker)
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

    def solve_dirichlet(self, f, h):
        """
        Solve with Dirichlet condition on deforming boundary
        """
        self.set_volume_forces(f)
        self.set_boundary_stress(h, self.deform_marker)
        bc = DirichletBC(self.V, h, self.mf, self.deform_marker)
        a = inner(self.sigma, grad(self.v))*dx
        L = inner(self.f,self.v)*dx
        solve(a==L, self.u_, bcs=bc)
        plot(self.u_)
        
if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("meshes/multimesh_1.xdmf") as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf_1.xdmf") as infile:
        infile.read(mvc, "name_to_read")
        mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    from create_meshes import inner_marker, outer_marker
    e_solve = ElasticitySolver(mesh, mf, free_marker=outer_marker,
                               deform_marker=inner_marker,
                               constant_mu=False)
    h = Expression(("x[0]","x[1]"), degree=1)
    f = Constant(("0","0"))
    e_solve.solve(f,h)

