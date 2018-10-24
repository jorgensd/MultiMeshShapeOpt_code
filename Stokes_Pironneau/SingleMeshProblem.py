from dolfin import (Mesh, MeshFunction, VectorElement, FiniteElement,
                    Function, FunctionSpace, XDMFFile, cpp, File,
                    MeshValueCollection, Constant, DirichletBC,
                    plot, TestFunctions, TrialFunctions,
                    inner, grad, dx, div, solve,
                    VectorFunctionSpace, BoundaryMesh, project,
                    Expression, dof_to_vertex_map,
                    vertex_to_dof_map, Vertex, Identity, tr,
                    Measure, TrialFunction, TestFunction, sym, ALE,
                    MeshEntity, SpatialCoordinate, as_vector, assemble)
from IPython import embed
import numpy
from pdb import set_trace
from matplotlib.pyplot import show

class StokesSolver():
    def __init__(self, mesh, mf, mesh_b, bc_dict, move_dict):
        self.mesh = mesh
        self.mesh_b = mesh_b
        self.backup = mesh.coordinates()
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
        self._init_geometric_functions()
        self.vfac = 1e4
        self.bfac = 1e2
        
    def _init_geometric_functions(self):
        x = SpatialCoordinate(self.mesh)
        VolOmega = assemble(Constant(1)*dx(domain=self.mesh))
        self.Vol0 = Constant(1 - VolOmega)
        self.bx0 = Constant((Constant(1./2)-assemble(x[0]*dx))/self.Vol0)
        self.by0 = Constant((Constant(1./2)-assemble(x[1]*dx))/self.Vol0)

    def compute_volume_bary(self):
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
        Solve Stokes problem
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
        #mesh.coordinates = self.backup
        #volume_function = self.InjectFunctionFromSurface(surface_deformation)
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

        # If we ever use mesh nodes as direct control
        # dStress = Measure("ds", subdomain_data=self.mf)
        # for marker in self.move_dict["Deform"]:
        #     L += inner(volume_function, v)*dStress(marker)
        
        bcs = []
        for marker in self.move_dict["Fixed"]:
            bcs.append(DirichletBC(self.S,
                                   Constant([0]*mesh.geometric_dimension()),
                                   self.mf, marker))
        s = Function(self.S)
        solve(a==L, s, bcs=bcs)
        self.perturbation = s


    def steepest_descent_update(self, step):
        s_descent = self.perturbation.copy(deepcopy=True)
        s_descent.vector()[:] *= step
        ALE.move(self.mesh, s_descent)

    def eval_J(self):
        # Evaluates J with current mesh
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
        print(step, J_eps-self.J, self.dJ)
        return J_eps

    # Helper function for moola
    def phi_dphi0(self):
        """ Return J and dJ with current mesh"""
        return self.J, self.dJ

    def eval_dJ(self):
        self.generate_mesh_deformation()
        self.dJ = assemble(self.dJ_form).inner(self.perturbation.vector())
        return self.dJ

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
        

    def InjectFunctionFromSurface(self, f):
        """ Take a CG1 function f defined on a surface mesh and return a 
        volume vector with same values on boundary but zero in volume
        """
        MaxDim = self.mesh.geometric_dimension()
        SpaceS = FunctionSpace(mesh, "CG", 1)
        SpaceV = VectorFunctionSpace(mesh, "CG", 1, MaxDim)
        F = Function(SpaceV)
        LocValues = numpy.zeros(F.vector().local_size())
        map = self.mesh_b.entity_map(0)
        OwnerRange = SpaceV.dofmap().ownership_range()
        d2v = dof_to_vertex_map(FunctionSpace(self.mesh_b, "CG", 1))
        v2d = vertex_to_dof_map(SpaceS)
        for i in range(int(f.vector().local_size()/MaxDim)):

            GVertID = Vertex(self.mesh_b, d2v[i]).index()
            PVertID = map[GVertID]
            PDof = v2d[PVertID]
            l_2_g_index = SpaceV.dofmap().local_to_global_index(PDof)
            IsOwned = (OwnerRange[0]/MaxDim <= l_2_g_index
                       and l_2_g_index<=OwnerRange[1]/MaxDim)
            if IsOwned:
                for j in range(MaxDim):
                    value = f.vector()[MaxDim*i+j]
                    LocValues[PDof*MaxDim+j] = value
        F.vector().set_local(LocValues)
        F.vector().apply("")
        return F


    def ReduceFunctionToSurface(self, volume_function):
        """
        Reduces a CG-1 function from a mesh to a CG-1 function on the boundary 
        mesh
        """
        
        MaxDim = self.mesh.geometric_dimension()
        surface_space = VectorFunctionSpace(self.mesh_b, "CG", 1)
        surfacevector = Function(surface_space)
        (sdmin, sdmax) = surfacevector.vector().local_range()
        sdmin = int(sdmin/MaxDim)
        sdmax = int(sdmax/MaxDim)
        LocValues = numpy.zeros(MaxDim*(sdmax-sdmin))
        
        VGlobal = numpy.zeros(len(volume_function.vector()))
        (vdmin, vdmax) = volume_function.vector().local_range()
        vdmin = int(vdmin/MaxDim)
        vdmax = int(vdmax/MaxDim)
        Own_min, Own_max = volume_function.function_space().dofmap().ownership_range()
        DofToVert = dof_to_vertex_map(FunctionSpace(self.mesh, "CG", 1))

        for i in range(vdmax-vdmin):
            Vert = MeshEntity(self.mesh, 0, DofToVert[i])
            GlobalIndex = Vert.global_index()
            for j in range(MaxDim):
                value = volume_function.vector()[MaxDim*i+j]
                VGlobal[MaxDim*GlobalIndex+j] = value
        mapa = self.mesh_b.entity_map(0)
        OwnerRange = surface_space.dofmap().ownership_range()
        DofToVert = dof_to_vertex_map(FunctionSpace(self.mesh_b, "CG", 1))
        for i in range(sdmax-sdmin):
            VolVert = MeshEntity(self.mesh, 0, mapa[int(DofToVert[i])])
            GlobalIndex = VolVert.global_index()

            for j in range(MaxDim):
                value = VGlobal[MaxDim*GlobalIndex+j]
                LocValues[MaxDim*i+j] = value

        surfacevector.vector().set_local(LocValues)
        surfacevector.vector().apply('')
        File("output/surface_movement.pvd") << surfacevector
        return surfacevector

        
if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("meshes/singlemesh.xdmf") as infile:
        infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    mesh_b = BoundaryMesh(mesh, "exterior")
    # S_b = VectorFunctionSpace(mesh_b, "CG", 1)
    # s_b = Function(S_b)
    #project(Expression(("cos(x[0])", "x[1]"), degree=1), S_b)

    from create_meshes import inflow, outflow, walls
    solver = StokesSolver(mesh, mf, mesh_b,
                          {inflow: Constant((1.0,0.0)),
                           walls: Constant((0.0,0.0))},
                          {"Fixed": [inflow, outflow], "Deform": [walls]})
    Js = []
    dJs = []
    
    def steepest_descent():
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
            meshfile << solver.mesh
            Ji = solver.eval_J()
            dJi = solver.eval_dJ()
            Js.append(Ji)
            dJs.append(dJi)
            if i>0:
                rel_red = (abs(Js[-1]-Js[-2])/Js[-1])
                print("Rel reduction: %.2e" % (rel_red))
                if rel_red < rel_tol:
                    break
            line_step = linesearch.search(solver.phi, None, solver.phi_dphi0())
            print("Step %.2e, J %.2e dJ %.2e" % (line_step, Ji, dJi))
            solver.steepest_descent_update(line_step)
            solver.generate_mesh_deformation()

    steepest_descent()
