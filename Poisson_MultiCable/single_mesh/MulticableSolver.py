from dolfin import *
from create_mesh import *
from IPython import embed
import numpy
class MultiCable():
    def __init__(self,positions, lmb_core,lmb_iso, lmb_fill, fs,res=0.01, state=None):
        if state is None:
            self.state = File("output/state.pvd")
        else:
            self.state = state
        self.J = 0
        self.dJ = 0
        self.num_cables = int(len(positions)/2)
        self.outT = File("output/T.pvd")
        self.lmb_core = lmb_core
        self.lmb_iso = lmb_iso
        self.lmb_fill = lmb_fill
        self.fs = fs
        self.init_mesh(positions,res=res)
        self.init_source_and_heat_coeff(fs, lmb_core, lmb_iso, lmb_fill)
        self.T = Function(self.V, name="Temperature")
        self.adjT = Function(self.V, name="Adjoint")
        self.q = 3 # Power of Functional
        self.obj = (1./self.q)*pow(abs(self.T), self.q)*dx # Functional
        self.objdT = self.T*pow(abs(self.T), self.q-2)
        self.T_amb = Constant(3.2) # Ambient Temperature
        self.c = Constant(0.01) # Reaction coefficient
        
    def update_mesh(self,positions):
        self.__init__(positions, self.lmb_core, self.lmb_iso,
                      self.lmb_fill, self.fs, self.res, self.state)
    
    def init_mesh(self, positions,res):
        self.res = res
        create_multicable(positions,res)
        self.mesh = Mesh()
        with XDMFFile("multicable.xdmf") as infile:
            infile.read(self.mesh)
        mvc = MeshValueCollection("size_t", self.mesh, 1)
        with XDMFFile("mf.xdmf") as infile:
            infile.read(mvc, "name_to_read")
        self.mf = cpp.mesh.MeshFunctionSizet(self.mesh, mvc)
        mvc = MeshValueCollection("size_t", self.mesh, 2)
        with XDMFFile("cf.xdmf") as infile:
            infile.read(mvc, "name_to_read")
        self.cf = cpp.mesh.MeshFunctionSizet(self.mesh, mvc)


    def init_source_and_heat_coeff(self, sources, metal, iso, fill):
        """ Initialize heat coefficient and source function 
        for all cables
        """
        if isinstance(metal,float) and isinstance(iso,float):
            # Wrapping single parameters as list
            metal = metal*numpy.ones(len(sources))
            iso = iso*numpy.ones(len(sources))
        self.V = FunctionSpace(self.mesh, "CG", 1)
        W = FunctionSpace(self.mesh, "DG", 0)
        self.f = Function(W)
        for i in range(self.num_cables):
            self.f.vector()[:] = (rubber_marker < self.cf.array())*(sources[i])
        self.lmb = Function(W)
        for i in range(self.num_cables):
            self.lmb.vector()[:] = ((self.cf.array() == rubber_marker)*iso[i] + 
                                    (self.cf.array() == fill_marker)*fill +
                                    (self.cf.array() == metal_marker)*metal[i])
    def alpha_heat_transfer(self, T):
        return Constant(1.0)

    def WeakCableShapeGradSurf(self, T, adjT, lmb, c, f, n):
        """ The Riez representer of the shape-surface gradient at an interface
        """
        # ("-") are exterior quantitites, [ ]_{+-}=-jump( )
        grad_tau_adjT = grad(adjT("-"))-dot(grad(adjT("-")),n)*n
        grad_tau_T = grad(T("-"))-dot(grad(T("-")),n)*n
        dJ = inner(grad_tau_T, grad_tau_adjT)*jump(lmb)\
             - adjT("-")*(jump(f)
             )\
             - lmb("-")*dot(n,grad(adjT("-")))*dot(jump(grad(T)),n)

        return -dJ

    def eval_J(self, positions):
        self.update_mesh(positions)
        # Evaluate J at current position
        n = FacetNormal(self.mesh)
        v = TestFunction(self.V)
        Ttmp = TrialFunction(self.V)
        constraint = inner(self.lmb*grad(Ttmp), grad(v))*dx \
                     -self.f*v*dx -self.c*v*Ttmp*dx
        constraint += self.alpha_heat_transfer(Ttmp)*(Ttmp-self.T_amb)*v*ds

        with Timer("USER_TIMING: Assemble State") as t:
            A = assemble(lhs(constraint))
            b = assemble(rhs(constraint))
        self.T.vector()[:]=0
        with Timer("USER_TIMING: Solve State") as t:
            solve(A, self.T.vector(), b, 'lu')
        self.J = assemble(self.obj)
        return self.J

    def save_state(self):
        self.state << self.T
    
    # Evaluate the shape gradient
    def eval_dJ(self, positions):
        self.eval_J(positions)
        dJ = []
        # Solve adjoint equation
        adj = TrialFunction(self.V)
        v = TestFunction(self.V)
        n = FacetNormal(self.mesh)
        constraint = inner(self.lmb*grad(adj), grad(v))*dx -self.c*v*adj*dx
        # FIXME: Add derivative of alpha in ext bc constraint,only works for alpha=1
        constraint += adj*1*v*ds# alpha_heat_transfer(T)*v*ds
        constraint += self.objdT*v*dx
        A = assemble(lhs(constraint))
        b = assemble(rhs(constraint))
        solve(A, self.adjT.vector(), b, 'lu')


        # Alternative to facet normal from femorph
        # from femorph import VolumeNormal
        # normal = VolumeNormal(self.mesh, [0], self.mf, [metaliso,isofill])
        # File("n.pvd") << normalu
        # normal = FacetNormal(cable_mesh)("-") # Outwards pointing normal
        dJ_Surf = self.WeakCableShapeGradSurf(self.T, self.adjT,
                                              self.lmb, self.c, self.f,
                                              n=n("-"))
        for i in range(self.num_cables):
            dSc1 = Measure("dS", subdomain_data=self.mf, subdomain_id=metaliso+i)
            dSc2 = Measure("dS", subdomain_data=self.mf, subdomain_id=isofill+i)
            # note volume normal is in opposite direction of MultiMesh
            gradx = assemble(-n[0]("+")*dJ_Surf*dSc1
                             -n[0]("+")*dJ_Surf*dSc2
                             + Constant(0)*dx(domain=self.mesh,
                                              subdomain_data=self.cf))
            grady = assemble(-n[1]("+")*dJ_Surf*dSc1
                             -n[1]("+")*dJ_Surf*dSc2
                             + Constant(0)*dx(domain=self.mesh,
                                              subdomain_data=self.cf))
            dJ.append(gradx)
            dJ.append(grady)
        self.dJ = numpy.array(dJ)
        return self.dJ

        
def convergence_rates(E_values, eps_values):
    r = []
    for i in range(1, len(eps_values)):
        r.append(numpy.log(E_values[i]/E_values[i-1])/
                 numpy.log(eps_values[i]/eps_values[i-1]))
    return r

def compute_angles(cable_positions):
    """
    Compute angles between three cables
    """
    c1 = cable_positions[0:2]
    c2 = cable_positions[2:4]
    c3 = cable_positions[4:6]
    a1 = (numpy.arccos(numpy.dot(c1-c2,c1-c3)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c1-c3,c1-c3))))/(2*numpy.pi)*360)
    a2 = (numpy.arccos(numpy.dot(c1-c2,c3-c2)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)
    a3 = (numpy.arccos(numpy.dot(c3-c1,c3-c2)/
                      (numpy.sqrt(numpy.dot(c3-c1,c3-c1)
                                  *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)
    print("Angles between the three cables: %.2f, %.2f, %.2f" %(a1,a2,a3))

if __name__ == "__main__":
    lmb_metal = 205.   # Heat coefficient aluminium
    lmb_iso = 0.03 # Heat coefficient of plastic
    lmb_air = 0.33   # Heat coefficient of brick
    c1 = numpy.array([0, 0.45])
    c2 = numpy.array([-0.4, -0.15])
    c3 = numpy.array([0.2,-0.4])
    cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
    compute_angles(cable_positions)
    sources = numpy.array([10,10,10])
    MC = MultiCable(cable_positions,lmb_metal,lmb_iso,lmb_air,sources,0.0176)
    tmp_out = File("output/mesh.pvd")
    J = MC.eval_J(cable_positions)
    print("----Number of cells---")
    print(MC.mesh.num_cells())
    tmp_out << MC.mesh
    MC.save_state()
    print(J)
    from Optimization import MultiCableOptimization
    # opt = MultiCableOptimization(int(len(cable_positions)/2), MC.eval_J, MC.eval_dJ)
    # sol = opt.solve(cable_positions)
    # compute_angles(sol)
    # MC.eval_J(sol)
    # tmp_out << MC.mesh
    list_timings(TimingClear.keep, [TimingType.wall])

