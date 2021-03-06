from dolfin import (Mesh, MeshFunction, VectorElement, FiniteElement,
                    Function, FunctionSpace, XDMFFile, cpp, File,
                    MeshValueCollection, Constant, DirichletBC,
                    plot, TestFunctions, TrialFunctions,
                    inner, grad, dx, div, solve,
                    VectorFunctionSpace, BoundaryMesh, project,
                    Expression, dof_to_vertex_map,
                    vertex_to_dof_map, Vertex, Identity, tr,
                    Measure, TrialFunction, TestFunction, sym, ALE,
                    MeshEntity, SpatialCoordinate, as_vector, assemble,
                    set_log_level, LogLevel)
set_log_level(LogLevel.ERROR)
from IPython import embed
import numpy
from pdb import set_trace
from matplotlib.pyplot import show

sm_def = File("output/sm_deform.pvd")
class StokesSolver():
    def __init__(self, mesh, mf, bc_dict, move_dict):
        """
        Inititalize Stokes solver with a mesh, its corresponding facet function,
        a dictionary describing boundary conditions and 
        a dictionary describing which boundaries are fixed in the shape optimization setting
        """
        self.mesh = mesh
        self.backup = mesh.coordinates().copy()
        self.mf = mf
        V2 = VectorElement("CG", mesh.ufl_cell(), 2)
        S1 = FiniteElement("CG", mesh.ufl_cell(), 1)
        TH = V2 * S1
        self.VQ = FunctionSpace(self.mesh, TH)
        self.w = Function(self.VQ)
        self.f = Constant([0.]*mesh.geometric_dimension())
        self.S = VectorFunctionSpace(self.mesh, "CG", 1)
        self.move_dict = move_dict
        self.J = 0
        self._init_bcs(bc_dict)
        self.solve()
        self.vfac = 1e4
        self.bfac = 1e2
        self._init_geometric_functions()
        self.eval_current_J()
        self.eval_current_dJ()
        self.outfile = File("output/u_singlemesh.pvd")
        self.outfile << self.u
        self.gradient_scale = 1
  
        self.iteration_counter = 1
        self.create_mapping_for_moving_boundary()
        
    def _init_geometric_functions(self):
        """
        Helper initializer to compute original volume and barycenter
        of the obstacle.
        """
        x = SpatialCoordinate(self.mesh)
        VolOmega = assemble(Constant(1)*dx(domain=self.mesh))
        self.Vol0 = Constant(1 - VolOmega)
        self.bx0 = Constant((Constant(1./2)-assemble(x[0]*dx))/self.Vol0)
        self.by0 = Constant((Constant(1./2)-assemble(x[1]*dx))/self.Vol0)

    def compute_volume_bary(self):
        """
        Compute barycenter and volume of obstacle with current mesh
        """
        x = SpatialCoordinate(self.mesh)
        VolOmega = assemble(Constant(1)*dx(domain=self.mesh))
        self.Vol = Constant(1 - VolOmega)
        self.bx = Constant((Constant(1./2)-assemble(x[0]*dx))/self.Vol)
        self.by = Constant((Constant(1./2)-assemble(x[1]*dx))/self.Vol)
        self.bx_off = self.bx - self.bx0
        self.by_off = self.by - self.by0
        self.vol_off = self.Vol - self.Vol0

    def _init_bcs(self, bc_dict):
        """
        Initialize boundary conditions for inlets and walls
        """
        self.bcs = []
        for marker in bc_dict.keys():
            bc = DirichletBC(self.VQ.sub(0), bc_dict[marker], self.mf, marker)
            self.bcs.append(bc)

    def solve(self):
        """
        Solve Stokes problem with current mesh
        """
        (u, p) = TrialFunctions(self.VQ)
        (v, q) = TestFunctions(self.VQ)
        a = inner(grad(u), grad(v))*dx - div(u)*q*dx - div(v)*p*dx
        l = inner(self.f, v)*dx
        solve(a==l, self.w, bcs=self.bcs)
        self.u, self.p = self.w.split()

        
    def generate_mesh_deformation(self):
        """
        Generates an linear elastic mesh deformation using the steepest
        gradient as stress on the boundary.
        """
        u, v = TrialFunction(self.S), TestFunction(self.S)

        def compute_mu(constant=True):
            """
            Compute mu as according to arxiv paper
            https://arxiv.org/pdf/1509.08601.pdf
            """
            mu_min=Constant(1)
            mu_max=Constant(500)
            if constant:
                return mu_max
            else:
                V = FunctionSpace(self.mesh, "CG",1)
                u, v = TrialFunction(V), TestFunction(V)
                a = inner(grad(u),grad(v))*dx
                l = Constant(0)*v*dx
                bcs = []
                for marker in self.move_dict["Fixed"]:
                    bcs.append(DirichletBC(V, mu_min, self.mf, marker))
                for marker in self.move_dict["Deform"]:
                    bcs.append(DirichletBC(V, mu_max, self.mf, marker))
                mu = Function(V)
                solve(a==l, mu, bcs=bcs)
                return mu
        mu = compute_mu(False)
        
        def epsilon(u):
            return sym(grad(u))
        def sigma(u,mu=500, lmb=0):
            return 2*mu*epsilon(u) + lmb*tr(epsilon(u))*Identity(2)
        a = inner(sigma(u,mu=mu), grad(v))*dx
        L = inner(Constant((0,0)), v)*dx
        L -= self.dJ_form
        
        bcs = []
        for marker in self.move_dict["Fixed"]:
            bcs.append(DirichletBC(self.S,
                                   Constant([0]*mesh.geometric_dimension()),
                                   self.mf, marker))
        s = Function(self.S)
        solve(a==L, s, bcs=bcs)
        self.perturbation = s


    def steepest_descent_update(self, step, out=False):
        """
        Updates the mesh in the steepest descent direction with steplength 'step'
        """
        s_descent = self.perturbation.copy(deepcopy=True)
        if out:
            sm_def << s_descent
        s_descent.vector()[:] *= step
        ALE.move(self.mesh, s_descent)

    def update_mesh_from_boundary_nodes(self, perturbation):
        """
        Deform mesh with boundary perturbation (a numpy array)
        given as Neumann input in an linear elastic mesh deformation.
        """
        # Reset mesh
        self.mesh.coordinates()[:] = self.backup
        volume_function = Function(self.S)

        for i in self.design_map.keys():
            volume_function.vector()[self.design_map[i]] = perturbation[i]
        u, v = TrialFunction(self.S), TestFunction(self.S)

        def compute_mu(constant=True):
            """
            Compute mu as according to arxiv paper
            https://arxiv.org/pdf/1509.08601.pdf
            """
            mu_min=Constant(1)
            mu_max=Constant(500)
            if constant:
                return mu_max
            else:
                V = FunctionSpace(self.mesh, "CG",1)
                u, v = TrialFunction(V), TestFunction(V)
                a = inner(grad(u),grad(v))*dx
                l = Constant(0)*v*dx
                bcs = []
                for marker in self.move_dict["Fixed"]:
                    bcs.append(DirichletBC(V, mu_min, self.mf, marker))
                for marker in self.move_dict["Deform"]:
                    bcs.append(DirichletBC(V, mu_max, self.mf, marker))
                mu = Function(V)
                solve(a==l, mu, bcs=bcs)
                return mu
        mu = compute_mu(False)
        
        def epsilon(u):
            return sym(grad(u))
        def sigma(u,mu=500, lmb=0):
            return 2*mu*epsilon(u) + lmb*tr(epsilon(u))*Identity(2)
        a = inner(sigma(u,mu=mu), grad(v))*dx
        L = inner(Constant((0,0)), v)*dx

        
        bcs = []
        for marker in self.move_dict["Fixed"]:
            bcs.append(DirichletBC(self.S,
                                   Constant([0]*mesh.geometric_dimension()),
                                   self.mf, marker))
            
        dStress = Measure("ds", subdomain_data=self.mf)
        from femorph import VolumeNormal
        n = VolumeNormal(mesh)

        # Enforcing node movement through elastic stress computation
        # NOTE: This strategy does only work for the first part of the deformation
        for marker in self.move_dict["Deform"]:
            L += inner(volume_function, v)*dStress(marker)
        # Direct control of boundary nodes
        # for marker in self.move_dict["Deform"]:
        #     bcs.append(DirichletBC(self.S, volume_function, self.mf, marker))
        s = Function(self.S)
        solve(a==L, s, bcs=bcs)
        self.perturbation = s
        ALE.move(self.mesh, self.perturbation)

    def eval_scipy_J(self, perturbation):
        """
        Evaluate functional with mesh perturbed with perturbation as input
        """
        # Update_mesh
        self.update_mesh_from_boundary_nodes(perturbation)
        J_eps = self.eval_current_J()
        return J_eps

    def eval_scipy_dJ(self, perturbation):
        """
        Evaluate functional gradient with mesh perturbation as input
        """
        # Evaluate gradient with given perturbation
        # Update mesh with perturbation
        self.update_mesh_from_boundary_nodes(perturbation)
        # Evaluate state equation for gradient input
        self.eval_current_J()
        # Update Gradient form with perturbation given by the elasticity
        # equations
        self.dJ = assemble(self.dJ_form)
        tmp_dJ = Function(self.S)
        tmp_dJ.vector()[:] = self.dJ

        # Reduce the problem to a boundary problem
        for i in self.design_map.keys():
            self.dJ_array[i] = self.gradient_scale*tmp_dJ.vector().get_local()[self.design_map[i]]
        return self.dJ_array

    def update_mesh_coordinates(self, perturbation):
        self.update_mesh_from_boundary_nodes(perturbation)
        self.backup = self.mesh.coordinates().copy()
    
    def eval_current_J(self):
        """
        Evaluates J with current mesh
        """
        self.solve()
        self.compute_volume_bary()
        self.J = assemble(inner(grad(self.u), grad(self.u))*dx)
        self.J += float(self.vfac*self.vol_off**2)
        self.J += float(self.bfac*self.bx_off**2)
        self.J += float(self.bfac*self.by_off**2)
        # Update gradient expression
        self.recompute_dJ()
        return self.J

    # Helper function for moola
    def phi(self, step):
        """
        Evaluate functional in steepest descent direction
        """
        # Perturbes the mesh and evaluates at current point
        self.backup = self.mesh.coordinates().copy()
        self.steepest_descent_update(step)
        self.solve()
        self.compute_volume_bary()
        J_eps = assemble(inner(grad(self.u), grad(self.u))*dx)
        J_eps += float(self.vfac*self.vol_off**2)
        J_eps += float(self.bfac*self.bx_off**2)
        J_eps += float(self.bfac*self.by_off**2)
        self.mesh.coordinates()[:] = self.backup
        return J_eps

    # Helper function for moola
    def phi_dphi0(self):
        """ Return J and dJ with current mesh"""
        return self.J, self.dJ

    def eval_current_dJ(self):
        """
        Returns gradient with current mesh
        """
        self.generate_mesh_deformation()
        self.dJ = assemble(self.dJ_form).inner(self.perturbation.vector())
        return self.dJ

    # FIXME: This could be simplified if we only store the integrand.
    def recompute_dJ(self):
        """
        Create gradient expression for deformation algorithm
        """

        # Recalculate barycenter and volume of obstacle
        self.compute_volume_bary()

        # Integrand of gradient
        x = SpatialCoordinate(mesh)
        dJ_stokes = -inner(grad(self.u), grad(self.u))
        dJ_vol = - Constant(2*self.vfac)*(self.Vol-self.Vol0)
        dJ_bar = Constant(2*self.bfac)/self.Vol*((self.bx-x[0])*self.bx_off
                                                 + (self.by-x[1])*self.by_off)
        self.integrand = dJ_stokes + dJ_vol + dJ_bar
        dDeform = Measure("ds", subdomain_data=self.mf)
        from femorph import VolumeNormal
        n = VolumeNormal(mesh)
        s = TestFunction(self.S)
        dJ = 0
        # Sum up over all moving boundaries
        for marker in self.move_dict["Deform"]:
            dJ += inner(s,n)*self.integrand*dDeform(marker)
        self.dJ_form = dJ
        
    def create_mapping_for_moving_boundary(self):
        """ Create a map from the boundary mesh vector function to the array
        of design variables
        """
        tmp_bcs = []
        for marker in self.move_dict["Deform"]:
            tmp_bcs.append(DirichletBC(self.S, Constant((1,1)), self.mf, marker))
        s_tmp = Function(self.S)
        [bc.apply(s_tmp.vector()) for bc in tmp_bcs]
        arr = s_tmp.vector().get_local()
        design_to_vec = {}
        j = 0
        for i in range(len(s_tmp.vector().get_local())):
            if arr[i]>0:
                design_to_vec[j] = i
                j+=1
        self.design_map = design_to_vec
        self.dJ_array = numpy.zeros(len(self.design_map.keys()))

    
    def callback(self, perturbation):
        print("Iteration %d" %self.iteration_counter)
        self.iteration_counter +=1
        self.mesh.coordinates()[:] = self.backup
        self.eval_scipy_dJ(perturbation)
        self.eval_current_J()
        print("Current J", self.J)

        self.outfile << self.u

        
if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("meshes/singlemesh.xdmf") as infile:
        infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)


    def scipy_optimization():
        """
        Using scipy and its optimization algorithms to solve the optimization problem
        """
        
        from create_meshes import inflow, outflow, walls
        solver = StokesSolver(mesh, mf,
                              {inflow: Constant((1.0,0.0)),
                               walls: Constant((0.0,0.0))},
                              {"Fixed": [inflow, outflow], "Deform": [walls]})
        print("original J", solver.J)
        from scipy.optimize import minimize
        # Design variables
        s_b = numpy.zeros(len(solver.design_map.keys()))
        solver.gradient_scale = 1e-4
        result = minimize(solver.eval_scipy_J, s_b,
                          jac=solver.eval_scipy_dJ,
                          method='BFGS',
                          callback=solver.callback,
                          options={'disp':True, 'maxiter':25})
        print("Update geometry and restart algorithm")
        # Update mesh to new perturbed mesh
        solver.update_mesh_coordinates(result['x'])

    #scipy_optimization()


    def steepest_descent():
        """
        Sovle stokes optimization problem with a simple steepest descent algorithm
        using armijo linesearch for steplength.
        """
        from create_meshes import inflow, outflow, walls
        solver = StokesSolver(mesh, mf,
                              {inflow: Constant((1.0,0.0)),
                               walls: Constant((0.0,0.0))},
                              {"Fixed": [inflow, outflow], "Deform": [walls]})
        Js = []
        dJs = []
        from moola.linesearch import ArmijoLineSearch
        max_it = 100
        start_stp = 2
        max_stp = 2
        min_stp = 1e-10
        red_tol = 1e-5
        linesearch = ArmijoLineSearch(start_stp=start_stp,
                                      stpmax=max_stp,
                                      stpmin=min_stp)
        meshfile = File("output/singlemesh.pvd")
        for i in range(max_it):
            print("Iteration %d" %i)
            meshfile << solver.mesh
            Ji = solver.eval_current_J()
            dJi = solver.eval_current_dJ()
            Js.append(Ji)
            dJs.append(dJi)
            if i>0:
                rel_red = (abs(Js[-1]-Js[-2])/Js[-1])
                print("Rel reduction: %.2e" % (rel_red))
                if rel_red < red_tol:
                    break
            line_step = linesearch.search(solver.phi, None, solver.phi_dphi0())
            print("Step %.2e, J %.2e dJ %.2e" % (line_step, Ji, dJi))
            solver.steepest_descent_update(line_step,out=True)
            solver.generate_mesh_deformation()
            
    steepest_descent()
