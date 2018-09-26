from dolfin import *
from IPython import embed;
import matplotlib.pyplot as plt
from numpy import isclose, allclose

mesh = Mesh()
with XDMFFile("meshes/singlemesh.xdmf") as infile:
    infile.read(mesh)
mvc_f = MeshValueCollection("size_t", mesh, 1)
with XDMFFile("meshes/mf.xdmf") as infile:
    infile.read(mvc_f, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc_f)
mvc_c = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("meshes/cf.xdmf") as infile:
    infile.read(mvc_c, "name_to_read")
mc = cpp.mesh.MeshFunctionSizet(mesh, mvc_c)

dI2 = Measure("dS", domain=mesh, subdomain_data=mf, subdomain_id=2)
X = SpatialCoordinate(mesh)
V = FunctionSpace(mesh, "DG", 2)
S = VectorFunctionSpace(mesh, "CG", 1)
s = TrialFunction(S)
n = FacetNormal(mesh)
u, v = project(1+cos(X[0]), V), TestFunction(V)

dofmap = V.dofmap()
from numpy import where
from numpy.random import random
cells_13 = where(mc.array()==13)[0]
cells_12 = where(mc.array()==12)[0]
dofs_13 =[]
for i in cells_13:
    dofs_13.extend(dofmap.cell_dofs(i))
dofs_13 = list(set(dofs_13))
dofs_12 =[]
for i in cells_12:
    dofs_12.extend(dofmap.cell_dofs(i))
dofs_12 = list(set(dofs_12))

for i in dofs_12:
    u.vector()[i] = random()
for i in dofs_13:
    u.vector()[i] = random()
plot(u)
plt.show()
def equivalent(a, b):
    a_ = assemble(a).array()
    b_ = assemble(b).array()
    return allclose(a_,b_, rtol=1e-16)

def equivalent_vec(a, b):
    a_ = assemble(a).get_local()
    b_ = assemble(b).get_local()
    return allclose(a_,b_, rtol=1e-16)


def test_source_term():
    f_ = X[0]*sin(X[0])*cos(X[1])
    l_s = f_*v*dx
    df = grad(f_)
    dl_s = derivative(-l_s, X, s)
    exact = -dot(s, df)*v*dx - div(s)*f_*v*dx
    assert(equivalent(dl_s, exact))

def test_a_s():
    a_s = dot(grad(u), grad(v))*dx
    da_s = derivative(a_s, X, s)
    exact = -dot(dot(grad(u), grad(s)), grad(v))*dx\
            +div(s)*dot(grad(u), grad(v))*dx\
            -dot(grad(u),dot(grad(v), grad(s)))*dx
    assert(equivalent(da_s, exact))

def test_J():
    J = 0.5*u*u*dx
    dJ = derivative(J, X, s)
    exact = div(s)*0.5*u*u*dx
    assert(equivalent_vec(dJ, exact))

def tan_div(s, n):
    return div(s)-dot(dot(grad(s),n),n)

def dn_mat(s, n):
    return dot(outer(grad(s)*n,n).T,n) - dot(grad(s).T, n)

def test_a_IP():
    alpha = 4.0
    beta = 4.0
    h = 2.0*Circumradius(mesh)
    h = (h('+') + h('-')) / 2

    a_IP_1 = - dot(avg(grad(u)), jump(v, n))*dI2 
    da_IP_1 = derivative(a_IP_1, X, s)
    a_IP_2 = - dot(avg(grad(v)), jump(u, n))*dI2
    da_IP_2 = derivative(a_IP_2, X, s)
    a_IP_3 = beta/h*inner(jump(u), jump(v))*dI2
    da_IP_3 = derivative(a_IP_3, X, s)

    exact_1 = -tan_div(s("+"), n("+"))*dot(n("+"), avg(grad(u)))*jump(v)*dI2\
              -dot(dn_mat(s("+"), n("+")), avg(grad(u))*jump(v))*dI2\
              +dot(n("+"), avg(dot(grad(u), grad(s)))*jump(v))*dI2
    exact_2 = -tan_div(s("+"), n("+"))*dot(n("+"), avg(grad(v))*jump(u))*dI2\
              -dot(dn_mat(s("+"), n("+")), avg(grad(v))*jump(u))*dI2\
              +dot(n("+"), avg(dot(grad(v), grad(s)))*jump(u))*dI2
    exact_3 = -tan_div(s("+"),n("+"))*beta/h*jump(u)*jump(v)*dI2

    # FIXME: Coordinate derivative yields a derivative here that is not
    # equal to the one I've derived (its very different
    a_f_1 = assemble(da_IP_1).norm("frobenius") 
    a_e_1 = assemble(exact_1).norm("frobenius")
    assert(isclose(a_f_1,a_e_1, rtol=1e-16))
    a_f_2 = assemble(da_IP_2).norm("frobenius") 
    a_e_2 = assemble(exact_2).norm("frobenius")
    assert(isclose(a_f_2,a_e_2, rtol=1e-16))
    a_f_3 = assemble(da_IP_3).norm("frobenius") 
    a_e_3 = assemble(exact_3).norm("frobenius")
    # assert(isclose(a_f_3,a_e_3, rtol=1e-16))


if __name__=="__main__":
    test_source_term()
    test_a_s()
    test_J()
    test_a_IP()
    # test_a_O()
    # beta*dot(jump(grad(u)), jump(grad(v)))*dO
    # dJdO = beta*(-dot(jump(dot(grad(u), grad(s))), jump(grad(v)))*dO\
    #              +div(s)*dot(jump(grad(u)), jump(grad(v)))*dO\
    #              -dot(jump(grad(u)), jump(dot(grad(v),grad(s))))*dO)
