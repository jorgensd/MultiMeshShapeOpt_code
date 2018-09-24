from dolfin import *
from numpy import allclose

def test_coordiant_derivative_facetnormal():
    mesh = UnitSquareMesh(10,10)
    V = VectorFunctionSpace(mesh, "CG", 1)
    v = interpolate(Expression(("x[0]","sin(x[1])"), degree=2),V)
    n = FacetNormal(mesh)
    
    a = inner(n, v)*ds
    X = SpatialCoordinate(mesh)
    s = TestFunction(V)
    dadX = derivative(a, X,s)

    dadEx = assemble((div(s)-dot(dot(grad(s),n),n))*inner(n,v)*ds
                     +dot(dot(outer(grad(s)*n, n).T,n)- dot(grad(s).T, n),v)*ds)
    dA = assemble(dadX)
    assert(allclose(dA.get_local(), dadEx.get_local()))
