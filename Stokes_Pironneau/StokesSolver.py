from dolfin import (ALE, cpp,
                    dof_to_vertex_map, Vertex, vertex_to_dof_map,
                    Constant, DirichletBC, Expression, SpatialCoordinate,
                    Mesh, MeshValueCollection, MeshFunction,
                    MultiMesh, MultiMeshFunction, MultiMeshFunctionSpace,
                    MultiMeshDirichletBC, MultiMeshSubSpace,
                    XDMFFile, File,
                    Point,
                    VectorElement, FiniteElement,
                    plot,
                    TestFunctions, TrialFunctions,
                    dX, dI, dI,
                    inner, grad, div, solve,
                    set_log_level, LogLevel)
from IPython import embed
from pdb import set_trace
from matplotlib.pyplot import show
set_log_level(LogLevel.ERROR)
class StokesSolver():
    def __init__(self, meshes, facetfunctions, cover_points,
                 bc_dict, move_dict):
        """
        Solve the stokes problem with multiple meshes.
        Arguments:
           meshes: List of dolfin meshes, in the order they should be added 
                   to the multimesh
           facetfunctions: List of FacetFunctions corresponding to the
                   meshes above
           cover_points: 
                   dict where the key is the mesh that should get covered cells,
                   value is at which point auto_cover should start
           bc_dict: Dictionary describing boundary conditions for the Stokes
                    problem
           move_dict: Dictionary describing which node that will be fixed and 
                   which are design variables in the optimization problem
        """
        self.__init_multimesh(meshes, cover_points)
        self.mfs = facetfunctions
        self.move_dict = move_dict
        self.V2 = VectorElement("CG", meshes[0].ufl_cell(), 2)
        self.S1 = FiniteElement("CG", meshes[0].ufl_cell(), 1)
        self.VQ = MultiMeshFunctionSpace(self.multimesh, self.V2*self.S1)
        V = MultiMeshFunctionSpace(self.multimesh, self.V2)
        Q = MultiMeshFunctionSpace(self.multimesh, self.S1)
        self.__init_bcs(bc_dict)
        self.w = MultiMeshFunction(self.VQ, name="State")
        self.u = MultiMeshFunction(V, name="u")
        self.p = MultiMeshFunction(Q, name="p")

        self.f = Constant([0.]*self.multimesh.part(0).geometric_dimension())
        self.N = len(meshes)
        self.outu = [File("output/u_%d.pvd" %i) for i in range(self.N)]
        self.outp = [File("output/p_%d.pvd" %i) for i in range(self.N)]

        
        self.J = 0
        self.dJ = 0
        self.opt_it = 0

        
    def __init_multimesh(self, meshes, cover_points):
        multimesh = MultiMesh()
        self.backup = []
        for mesh in meshes:
            multimesh.add(mesh)
            self.backup.append(mesh.coordinates().copy())
        multimesh.build()
        for key in cover_points.keys():
            multimesh.auto_cover(key, cover_points[key])            
        self.multimesh = multimesh
        self.cover_points = cover_points

    def __init_bcs(self, bc_dict):
        """ 
        Initialize velocity dirichlet bcs on each mesh according to dictionary
        """
        V = MultiMeshSubSpace(self.VQ, 0)
        bcs = []
        for i in bc_dict:
            for marker in bc_dict[i]:
                bc = MultiMeshDirichletBC(V, bc_dict[i][marker],
                                          self.mfs[i], marker, i)
                bcs.append(bc)
        
if __name__ == "__main__":
    meshes = []
    for i in range(2):
        mesh_i = Mesh()
        with XDMFFile("meshes/multimesh_%d.xdmf" %i) as infile:
            infile.read(mesh_i)
        meshes.append(mesh_i)
    mfs = []
    for i in range(2):
        mvc = MeshValueCollection("size_t", meshes[i], 1)
        with XDMFFile("meshes/mf_%d.xdmf" %i) as infile:
            infile.read(mvc, "name_to_read")
        mfs.append(cpp.mesh.MeshFunctionSizet(meshes[i], mvc))

    cover = {0: Point(0.5,0.5)}
    from create_meshes import inflow, outflow, walls, inner_marker, outer_marker
    bc_dict = {0: {inflow: Constant((1.0,0.0)), walls: Constant((1,0))},
               1: {inner_marker: Constant((0,0))}}
    move_dict = {0: {"Fixed": [inflow, outflow, walls]},
                 1: {"Deform": [inner_marker, outer_marker]}}
    StokesSolver(meshes, mfs, cover, bc_dict, move_dict)
