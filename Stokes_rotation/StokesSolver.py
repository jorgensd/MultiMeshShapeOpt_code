from dolfin import *
from IPython import embed
from pdb import set_trace
import numpy as np
import matplotlib.pyplot as plt

class StokesSolver():

    def __init__(self, points, thetas, mesh_names, facet_func_names):
        """
        Initialize Stokes solver for objects located at "points"
        with orientation "thetas" with meshes from "mesh_names"
        and facet_functions from "facet_func_names".
        Arguments:
            points list(dolfin.Point) - The rotational centers
            thetas list(float)        - The initial orientations
            mesh_names list(str)      - The mesh-filenames
            facet_func_names list(str)- The facet_func filnames
        """

        self.N = len(points)
        self.points = points
        self.s = [Expression(("-x[1]+%s" % points[i][1],
                              "x[0]-%s" % points[i][0]), degree=3)
                  for i in range(self.N)]
        self.thetas = thetas
        self.init_multimesh(mesh_names, facet_func_names, self.thetas,
                            self.points)
        self.V2 = VectorElement("CG", triangle, 2)
        self.S1 = FiniteElement("CG", triangle, 1)
        self.VQ = MultiMeshFunctionSpace(self.multimesh, self.V2*self.S1)


    def init_multimesh(self, meshes_n, facet_funcs,theta,p): 
        multimesh = MultiMesh()
        mfs = []
        meshes = []
        for i in range(self.N+1):
            mesh_i = Mesh()
            with XDMFFile(meshes_n[i]) as infile:
                infile.read(mesh_i)
            mvc = MeshValueCollection("size_t", mesh_i, 1)
            with XDMFFile(facet_funcs[i]) as infile:
                infile.read(mvc, "name_to_read")
            mfs.append(cpp.mesh.MeshFunctionSizet(mesh_i, mvc))
            if i>0:
                mesh_i.translate(Point(p[i-1][0]-0.5, p[i-1][1]-0.5))
                mesh_i.rotate(theta[i-1], 2, p[i-1])
            meshes.append(mesh_i)
            multimesh.add(mesh_i)
        multimesh.build()
        self.mfs = mfs
        self.meshes = meshes
        self.multimesh = multimesh


points = [Point(0.5,0.25), Point(0.75,0.52)]
thetas = [90, 47]
pre = "meshes/"
meshes = [pre+"multimesh_0.xdmf", pre+"multimesh_1.xdmf",pre+"multimesh_1.xdmf"]
mfs = [pre+"mf_0.xdmf", pre+"mf_1.xdmf",pre+"mf_1.xdmf"]
ss = StokesSolver(points, thetas, meshes, mfs)
embed()
