from dolfin import *
from IPython import embed;
import matplotlib.pyplot as plt
from numpy import allclose

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

dInterface = Measure("dS", domain=mesh, subdomain_data=mf, subdomain_id=2)
X = SpatialCoordinate(mesh)
V = FunctionSpace(mesh, "CG", 1)
S = VectorFunctionSpace(mesh, "CG", 1)
s = TrialFunction(S)
n = FacetNormal(mesh)
u, v = project(1+cos(X[0]), V), TestFunction(V)



def equivalent(a, b):
    a_ = assemble(a).array()
    b_ = assemble(b).array()
    return allclose(a_,b_)

def equivalent_vec(a, b):
    a_ = assemble(a).get_local()
    b_ = assemble(b).get_local()
    return allclose(a_,b_)


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
    a_IP = - (dot(avg(grad(u)), jump(v, n))*dInterface 
              + dot(avg(grad(v)), jump(u, n))*dInterface)\
              +alpha/h*jump(u)*jump(v)*dInterface
    da_IP = derivative(a_IP, X, s)
    exact =-tan_div(s("+"), n("+"))*dot(n("+"),avg(grad(u)))*jump(v)*dInterface\
        -tan_div(s("+"), n("+"))*dot(n("+"), avg(grad(v))*jump(u))*dInterface\
        -tan_div(s("+"),n("+"))*beta/h*jump(u)*jump(v)*dInterface\
        -dot(dn_mat(s("+"), n("+")), avg(grad(v))*jump(u))*dInterface\
        -dot(dn_mat(s("+"), n("+")), avg(grad(u))*jump(v))*dInterface\
        +dot(n("+"), avg(dot(grad(u), grad(s)))*jump(v))*dInterface\
        +dot(n("+"), avg(dot(grad(v), grad(s)))*jump(u))*dInterface

    return equivalent(da_IP, exact)

if __name__=="__main__":
    test_source_term()
    test_a_s()
    test_J()
    test_a_IP()
# def a_O(u,v, beta=4.0):
#     return beta*dot(jump(grad(u)), jump(grad(v)))*dO


#     dJdO = -dot(jump(dot(grad(u), grad(s))), jump(grad(v)))*dO\
#            +div(s)*dot(jump(grad(u)), jump(grad(v)))*dO\
#            -dot(jump(grad(u)), jump(dot(grad(v),grad(s))))*dO
