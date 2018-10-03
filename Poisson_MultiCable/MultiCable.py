from dolfin import *
import numpy

class MultiCable():
    def __init__(self, scales, positions, lmb_core, lmb_iso, lmb_fill, fs):
        self.num_cables = len(scales)
        self.lmb_core = lmb_core
        self.lmb_iso = lmb_iso
        self.lmb_fill = lmb_fill
        self.init_multimesh(scales, positions)
        self.init_source_and_heat_coeff(fs, lmb_core, lmb_iso, lmb_fill)
        self.T = MultiMeshFunction(self.V, name="Temperature")
        self.adjT = MultiMeshFunction(self.V, name="Adjoint")
        self.alpha = 4.0 # First MultiMesh stab parameter
        self.beta = 4.0 # Second MultiMesh stab parameter
        self.q = 3 # Power of Functional
        self.obj = (1./self.q)*pow(abs(self.T), self.q)*dX # Functional
        self.objdT = self.T*pow(abs(self.T), self.q-2) # Derivative of functional integrand
        self.T_amb = Constant(3.2) # Ambient Temperature
        self.c = Constant(0.04) # Reaction coefficient

    def alpha_heat_transfer(self, T):
        return Constant(1.0)
        
    def init_multimesh(self,scales, positions):
        """
        Initializes a Multimesh with the scales of each interal cable and 
        initial position 
        """
        multimesh = MultiMesh()
        self.cable_subdomains = []
        self.cable_facets = []

        # Background mesh
        mesh_0 = Mesh("meshes/cable.xml")
        facets_0 = MeshFunction("size_t", mesh_0,
                                "meshes/cable_facet_region.xml")
        # Move origo-centered meshes to initial posiion
        self.cable_meshes = [mesh_0]
        self.num_cables = len(scales)
        for i in range(self.num_cables):
            cable = Mesh("meshes/inner_cable_halo.xml")
            cable.coordinates()[:,:]*= scales[i]
            cable.translate(Point(positions[2*i], positions[2*i+1]))
            self.cable_meshes.append(cable)
            subdomains = MeshFunction("size_t", cable,
                                      "meshes/inner_cable_halo_physical_region.xml")
            self.cable_subdomains.append(subdomains)
            facets = MeshFunction("size_t", cable,
                                  "meshes/inner_cable_halo_facet_region.xml")
            self.cable_facets.append(facets)

        # Create multimesh
        self.multimesh = MultiMesh()
        for cable in self.cable_meshes:
            self.multimesh.add(cable)
        self.multimesh.build()
        
    def init_source_and_heat_coeff(self, sources, metal, iso, fill):
        """ Initialize heat coefficient and source function 
        for all cables
        """
        if isinstance(metal,float) and isinstance(iso,float):
            # Wrapping single parameters as list
            metal = metal*numpy.ones(len(sources))
            iso = iso*numpy.ones(len(sources))
        self.V = MultiMeshFunctionSpace(self.multimesh, "CG", 1)
        W = MultiMeshFunctionSpace(self.multimesh, "DG", 0)
        W0 = FunctionSpace(self.cable_meshes[0], "DG", 0)
        self.f = MultiMeshFunction(W)
        self.f.assign_part(0, project(Constant(0.0), W0))
        for i in range(self.num_cables):
            X = FunctionSpace(self.cable_meshes[i+1], "DG", 0)
            fx = Function(X)
            fx.vector()[:] = (14.5 < self.cable_subdomains[i].array())*(sources[i])
            self.f.assign_part(i+1, fx)
        self.lmb = MultiMeshFunction(W)
        self.lmb.assign_part(0, project(Constant(fill), W0))
        plot(self.lmb.part(0))
        for i in range(self.num_cables):
            X = FunctionSpace(self.cable_meshes[i+1], "DG", 0)            
            lmbx = Function(X)
            lmbx.vector()[:] = ((self.cable_subdomains[i].array() == 13)*fill + 
                                (self.cable_subdomains[i].array() == 14)*iso[i] +
                                (self.cable_subdomains[i].array() == 15)*metal[i])
            self.lmb.assign_part(i+1, lmbx)
            
    def WeakCableShapeGradSurf(self, T, adjT, lmb, c, f, n):
        """ The Riez representer of the shape-surface gradient at an interface
        """
        # ("-") are exterior quantitites, [ ]_{+-}=-jump( )
        grad_tau_adjT = grad(adjT("-"))-dot(n,grad(adjT("-")))*n
        grad_tau_T = grad(T("-"))-dot(n,grad(T("-")))*n
        dJ = inner(grad_tau_T, grad_tau_adjT)*jump(lmb)\
             - adjT("-")*(jump(f)
             )\
             - lmb("-")*dot(n,grad(adjT("-")))*dot(n,jump(grad(T)))

        return -dJ

    
    def update_mesh(self,cable_positions):
        """ Translate all new_cables to a new center """
        cable_positions = cable_positions.reshape(-1, 2)
        for i, ((newx, newy), cable_mesh, cable_facet) in enumerate(
                zip(cable_positions, self.cable_meshes[1:], self.cable_facets)):
        # Move mesh to new positions
            dSc = Measure("dS", subdomain_data=cable_facet, subdomain_id=17)
            area = assemble(Constant(1)*dx(domain=cable_mesh))
            oldx = assemble(Expression("x[0]", degree=1)*dx(domain=cable_mesh))/area
            oldy = assemble(Expression("x[1]", degree=1)*dx(domain=cable_mesh))/area
            cable_mesh.translate(Point(-oldx+newx, -oldy+newy))
        self.multimesh.build()  # Rebuild the multimesh


    """ Evaluate the functional with given cable_positions"""
    def eval_J(self, cable_positions):
        # Update mesh
        self.update_mesh(cable_positions)
        n = FacetNormal(self.multimesh)
        h = 2.0*Circumradius(self.multimesh)
        h = (h('+') + h('-')) / 2
        v = TestFunction(self.V)
        Ttmp = TrialFunction(self.V)
        constraint = inner(self.lmb*grad(Ttmp), grad(v))*dX \
                     -self.f*v*dX -self.c*v*Ttmp*dX
        constraint += self.alpha_heat_transfer(Ttmp)*(Ttmp-self.T_amb)*v*ds
        constraint += - inner(avg(self.lmb*grad(Ttmp)), jump(v, n))*dI \
                      - inner(avg(self.lmb*grad(v)), jump(Ttmp, n))*dI \
                      + self.alpha/h*jump(Ttmp)*jump(v)*dI   \
                      + self.beta*self.lmb*inner(jump(grad(Ttmp)), jump(grad(v)))*dO
        A = assemble_multimesh(lhs(constraint))
        b = assemble_multimesh(rhs(constraint))
        self.T.vector()[:]=0
        self.V.lock_inactive_dofs(A, b)
        solve(A, self.T.vector(), b, 'lu')
    
        return assemble_multimesh(self.obj)


    # Evaluate the shape gradient
    def eval_dJ(self, cable_positions):
        # Update mesh
        self.eval_J(cable_positions)
        dJ = []
        # Solve adjoint equation
        adj = TrialFunction(self.V)
        v = TestFunction(self.V)
        n = FacetNormal(self.multimesh)
        h = 2.0*Circumradius(self.multimesh)
        h = (h('+') + h('-')) / 2
        constraint = inner(self.lmb*grad(adj), grad(v))*dX -self.c*v*adj*dX
        # FIXME: Add derivative of alpha in ext bc constraint,only works for alpha=1
        constraint += adj*1*v*ds# alpha_heat_transfer(T)*v*ds
        constraint += - inner(avg(self.lmb*grad(adj)), jump(v, n))*dI \
                      - inner(avg(self.lmb*grad(v)), jump(adj, n))*dI \
                      + self.alpha/h*jump(adj)*jump(v)*dI   \
                      + self.beta*self.lmb*inner(jump(grad(adj)),
                                                 jump(grad(v)))*dO
        constraint += self.objdT*v*dX
        A = assemble_multimesh(lhs(constraint))
        b = assemble_multimesh(rhs(constraint))
        self.V.lock_inactive_dofs(A, b)
        solve(A, self.adjT.vector(), b, 'lu')

        for i, (cable_mesh, cable_facet,cable_subdomain) in enumerate(
                zip(self.cable_meshes[1:], self.cable_facets,
                    self.cable_subdomains)):
            T_cable = self.T.part(i+1, deepcopy=True)
            adjT_cable = self.adjT.part(i+1, deepcopy=True)
            lmb_cable = self.lmb.part(i+1, deepcopy=True)
            f_cable = self.f.part(i+1, deepcopy=True)
            # Alternative to facet normal from femorph
            # from femorph import VolumeNormal
            # normal = VolumeNormal(cable_mesh, [0], cable_facet, [16,17])
            normal = FacetNormal(cable_mesh)("-") # Outwards pointing normal
            dJ_Surf = self.WeakCableShapeGradSurf(T_cable, adjT_cable,
                                                  lmb_cable, self.c, f_cable,
                                                  n=normal)
            dSc1 = Measure("dS", subdomain_data=cable_facet, subdomain_id=16)
            dSc2 = Measure("dS", subdomain_data=cable_facet, subdomain_id=17)
            gradx = assemble(normal[0]*dJ_Surf*dSc1
                             +normal[0]*dJ_Surf*dSc2
                             + Constant(0)*dx(domain=cable_mesh,
                                              subdomain_data=cable_subdomain))
            grady = assemble(normal[1]*dJ_Surf*dSc1
                             +normal[1]*dJ_Surf*dSc2
                             + Constant(0)*dx(domain=cable_mesh,
                                              subdomain_data=cable_subdomain))

            dJ.append(gradx)
            dJ.append(grady)
        return numpy.array(dJ)

    
if __name__ == "__main__":
    lmb_metal = 205.      # Heat coefficient aluminium
    lmb_insulation = 0.03 # Heat coefficient of plastic
    lmb_air = 0.15        # Heat coefficient of brick
    c1 = numpy.array([0, 0.45])
    c2 = numpy.array([-0.4, -0.15])
    c3 = numpy.array([0.2,-0.4])
    cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
    scales = numpy.array([1,1,1])   
    sources = numpy.array([5,5,5])
    
    MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                          lmb_air, sources)
    J_T = MC.eval_J(cable_positions)
    dJ_T = MC.eval_dJ(cable_positions)
    print(J_T)
    print(dJ_T)
