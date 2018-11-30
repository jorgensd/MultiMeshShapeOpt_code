from dolfin import *
from dolfin_adjoint import *
from create_mesh import *
import matplotlib.pyplot as plt
from IPython import embed
import numpy
set_log_level(LogLevel.CRITICAL)

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
        self.s = Function(VectorFunctionSpace(self.mesh, "CG", 1))
        self.init_source_and_heat_coeff(fs, lmb_core, lmb_iso, lmb_fill)
        self.update_mesh()
        self.T = Function(self.V, name="Temperature")
        self.adjT = Function(self.V, name="Adjoint")
        self.q = 3 # Power of Functional
        self.obj = (1./self.q)*pow(abs(self.T), self.q)*dx
        self.objdT = self.T*pow(abs(self.T), self.q-2)
        self.T_amb = Constant(3.2) # Ambient Temperature
        self.c = Constant(0.01) # Reaction coefficient

    def update_mesh(self):
        ALE.move(self.mesh, self.s)

        
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

    def eval_J(self):
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


if __name__ == "__main__":
    lmb_metal = 205.   # Heat coefficient aluminium
    lmb_iso = 0.03 # Heat coefficient of plastic
    lmb_air = 0.33   # Heat coefficient of brick
    c1 = numpy.array([0, 0.45])
    cable_positions = numpy.array([c1[0],c1[1]])
    sources = numpy.array([10])
    MC = MultiCable(cable_positions,lmb_metal,lmb_iso,lmb_air,sources,0.0176)
    
    J = MC.eval_J()
    Jhat = ReducedFunctional(J, Control(MC.s))
    Jhat.optimize()
    tape.visualise("output/dot.dot", dot=True)
    s_tmp = Function(MC.s.function_space())
    for i in range(MC.num_cables):
        bc = DirichletBC(MC.s.function_space(), Constant((0,0.1)), MC.mf,
                         metaliso+i)
        bc.apply(s_tmp.vector())
        bc = DirichletBC(MC.s.function_space(), Constant((0,0.1)), MC.mf,
                         isofill+i)
        bc.apply(s_tmp.vector())            
    s_0 = Function(MC.s.function_space())

    taylor_test(Jhat, s_0, s_tmp)

