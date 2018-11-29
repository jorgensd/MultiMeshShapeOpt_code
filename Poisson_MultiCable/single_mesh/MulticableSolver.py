
from dolfin import *
from create_mesh import *
from IPython import embed


x0,y0 = 0,0
create_multicable(x0,y0)

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


embed()
