from dolfin import (Mesh, MeshFunction, VectorElement, FiniteElement,
                    Function, FunctionSpace, XDMFFile, cpp,
                    MeshValueCollection, Constant, DirichletBC)
from IPython import embed
from pdb import set_trace

class StokesSolver():
    def __init__(self, mesh, mf, bc_dict):
        self.mesh = mesh
        self.mf = mf
        V2 = VectorElement("CG", mesh.ufl_cell(), 2)
        S1 = FiniteElement("CG", mesh.ufl_cell(), 1)
        TH = V2 * S1
        self.VQ = FunctionSpace(self.mesh, TH)
        self.w = Function(self.VQ)
        self._init_bcs(bc_dict)
        self.f = Constant([0.]*mesh.geometric_dimension())
    def _init_bcs(self, bc_dict):
        """
        Initialize boundary conditions for inlets and walls
        """
        self.bcs = []
        for marker in bc_dict.keys():
            bc = DirichletBC(self.VQ.sub(0), bc_dict[marker], self.mf, marker)
            self.bcs.append(bc)

    def solve(self):
        (u, p) = TrialFunctions(VQ)
        (v, q) = TestFunctions(VQ)
        

if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("meshes/singlemesh.xdmf") as infile:
            infile.read(mesh)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("meshes/mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    from create_meshes import inflow, outflow, walls
    StokesSolver(mesh, mf, {inflow: Constant((1.0,0.0)),
                            walls: Constant((0.0,0.0))})
