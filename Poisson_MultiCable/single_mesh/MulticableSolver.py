from dolfin import *
from create_mesh import *
from IPython import embed



class MultiCable():
    def __init__(self,positions, lmb_iso, lmb_core, lmb_fill, fs):
        self.J = 0
        self.dJ = 0
        self.num_cables = int(len(positions)/2)
        self.outT = File("output/T.pvd")
        self.lmb_core = lmb_core
        self.lmb_iso = lmb_iso
        self.lmb_fill = lmb_fill
        self.init_mesh(positions)
        # self.init_source_and_heat_coeff(fs, lmb_core, lmb_iso, lmb_fill)
        # self.T = MultiMeshFunction(self.V, name="Temperature")
        # self.adjT = MultiMeshFunction(self.V, name="Adjoint")
        # self.alpha = 4.0 # First MultiMesh stab parameter
        # self.beta = 4.0 # Second MultiMesh stab parameter
        # self.q = 3 # Power of Functional
        # self.obj = (1./self.q)*pow(abs(self.T), self.q)*dX # Functional
        # self.objdT = self.T*pow(abs(self.T), self.q-2) # Derivative of functional integrand
        # self.T_amb = Constant(3.2) # Ambient Temperature
        # self.c = Constant(0.01) # Reaction coefficient

    def init_mesh(self, positions):
        create_multicable(positions)
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

cable_positions = [0,0,0.5,0.5]
MC = MultiCable(cable_positions,0,0,0,0)

