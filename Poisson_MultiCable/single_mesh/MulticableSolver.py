
from dolfin import *
from create_mesh import *
from IPython import embed


cable_positions = [0,0,0.5,0.5]
create_multicable(cable_positions)

mesh = Mesh()
with XDMFFile("multicable.xdmf") as infile:
    infile.read(mesh)

mvc = MeshValueCollection("size_t", mesh, 1)
with XDMFFile("mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("cf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# class MultiCable():
#     (self, scales, positions, lmb_core, lmb_fill, fs):
    


embed()
