# Copyright (C) 2017 Jorgen Schartum Dokken
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2017-02-09
# Last changed: 2017-08-11
#
# This demo program solves the shape optimization problem of reduction of drag
# over a body subject to Stokes-flow. Overlapping meshes are used for simulation
# of the Stokes flow. Deformation of the mesh is done either with a Laplacian
# smoothing function or a projection equation. Central Voroi Tesselation smoothing
# is employed at the boundaries of the top mesh, to keep the mesh quality.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from dolfin import *
import moola
from IPython import embed
from pdb import set_trace
set_log_level(40)

"""
Stokes Solver for multimesh problem, with option of having solution saved as output.
   Input:
         multimesh - A multimesh consisting of two meshes, one background mesh
                     (mesh_0) with a hole where mesh_1 should be. mesh_1 should
                     have overlapping outer boundaries, while the inner boundary
                     is the flow obstacle.
         out       - A dictionary containing global the functionspace, solution
                     function and output files.
   Output:
         u,p       - MultiMeshFunctions containing velocity and pressure solutions of the stokes equation
"""
def StokesSolve(multimesh, out=None):
    if out==None:
        # Define Finite Element-spaces if this is a Armijo-linesearch.
        V2 = VectorElement("CG", triangle, 2)
        S1 = FiniteElement("CG", triangle, 1)
        TH = V2 * S1
        VQ = MultiMeshFunctionSpace(multimesh, TH)
        w = MultiMeshFunction(VQ)
    else:
        # Use Global spaces if this is a standard call
        VQ = out["VQ"]
        w = out["w"]

    [mf_0, mf_1] = load_facet_function(multimesh)
    
    # Define trial and test functions
    (u, p) = TrialFunctions(VQ)
    (v, q) = TestFunctions(VQ)
    f = Constant((0.0, 0.0))

    # Define facet normal and mesh size
    n = FacetNormal(multimesh)
    h = 2.0*Circumradius(multimesh)

    # Define bilinear and linear form
    a = a_h(u, v, n, h) + b_h(v, p, n) + b_h(u, q, n) + s_O(u, v) \
        + s_C(u, p, v, q, h)
    L  = l_h(v, q, f) + l_C(v, q, f, h)

    multimesh.build()
    multimesh.auto_cover(0, Point(0.5,0.5))

    # Create boundary conditions
    inflow_value = Expression(("1.0", "0.0"),degree=1)
    outflow_value = Constant(0)
    noslip_value = inflow_value
    obstacle_value = Expression(("0.0", "0.0"), degree=2)
    V = MultiMeshSubSpace(VQ, 0)
    Q = MultiMeshSubSpace(VQ, 1)
    #bc0 = MultiMeshDirichletBC(V, noslip_value,  mf_0, 3, 0)
    bc0 = MultiMeshDirichletBC(V,inflow_value,  mf_0, 3, 0)
    bc1 = MultiMeshDirichletBC(V, inflow_value,  mf_0, 1, 0)
    bc3 = MultiMeshDirichletBC(V, obstacle_value,  mf_1, 2 ,1)
    bcs = [bc0, bc1, bc3]

    # Assemble linear system, apply boundary conditions and solve
    assemble_time = -time.time()
    A = assemble_multimesh(a)
    b = assemble_multimesh(L)
    assemble_time += time.time()
    [bc.apply(A, b) for bc in bcs]
    VQ.lock_inactive_dofs(A, b)
    solve_time = -time.time()
    solve(A, w.vector(), b, "mumps")
    solve_time += time.time()
    print("Assemble time: %.2e" % assemble_time)
    print("Solve time: %.2e" % solve_time)
    # Splitting the mixed-multimeshfunction
    u, p = splitMMF(w, multimesh)
    return u, p


def Riez_Representation(mesh, step, direction, alpha=1):
    f=File("output/test2.pvd")
    W1 = VectorFunctionSpace(mesh, "CG", 1)
    w = TestFunction(W1)
    v = TrialFunction(W1)   
    n = FacetNormal(mesh)
    value = Expression("100*exp(-(x[0]-0.5)*(x[0]-0.5)/0.01)", degree=2)
    [mf_0, mf_1] = load_facet_function(multimesh)
    ds1 = ds(domain=mesh, subdomain_data=mf_1)
    alphas = [0, 1e-4, 1e-2, 1e-0, 1e2, 1e4, 1e6, 1e8]
    w1 = Function(W1)
    for alpha in alphas:
        from femorph import VolumeNormal
        normal = VolumeNormal(mesh, [0], facet)

        a = (Constant(alpha)*inner(grad(v),grad(w))+inner(v,w))*dx
        #l = step*inner(w, n*direction)*ds1(2)
        # l = step*inner(w, n*value)*ds1(2)
        l = inner(Constant((0.,0.)), v)*dx
        bc = DirichletBC(W1, normal*direction, facet, 2)
        solve(a ==l , w1,bcs=bc,solver_parameters={"linear_solver":"mumps"})
        f << w1

    return w1


def vector_norm(f):
    # computes the function sqrt(f[0]^2+f[1]^2)
    multimesh = f.function_space().multimesh()
    element = f.function_space().info[0]
    V_scalar = MultiMeshFunctionSpace(multimesh, element.family(),
                               element.degree())
    V_vector = f.function_space()
    f_x, f_y = MultiMeshFunction(V_scalar), MultiMeshFunction(V_scalar)
    for i in range(f.function_space().multimesh().num_parts()):
        f_i = f.part(i,deepcopy=True)
        fx_i, fy_i = f_x.part(i,deepcopy=True), f_y.part(i,deepcopy=True)
        f_i.vector()[:] *= f_i.vector()
        assigner_V_vector_to_V_scalar = FunctionAssigner([V_scalar.part(i),
                                                          V_scalar.part(i)],
                                                         V_vector.part(i))
        assigner_V_vector_to_V_scalar.assign([fx_i,fy_i], f_i)
        fx_i.vector().axpy(1, fy_i.vector())
        fx_i.vector().set_local(np.sqrt(fx_i.vector().get_local()))
        fx_i.vector().apply('')
        f_x.assign_part(i, fx_i)
    return f_x

def load_facet_function(multimesh):
    mfs = []
    for i in range(multimesh.num_parts()):
        mvc = MeshValueCollection("size_t", multimesh.part(i), 1)
        with XDMFFile("meshes/mf_%d.xdmf" %i) as infile:
            infile.read(mvc, "name_to_read")
        mfs.append(cpp.mesh.MeshFunctionSizet(multimesh.part(i), mvc))
    return mfs

def Laplacian(mesh, mf_1, n, step, direction, alpha=1e-2):
    # Smoothed H1 representation
    W1 = VectorFunctionSpace(mesh, "CG", 1)
    v = TestFunction(W1)
    w = TrialFunction(W1)
    a = inner(Constant(alpha)*grad(w),grad(v))*dx+inner(w,v)*dx
    l = inner(Constant((0,0)),v)*dx
    move = step*n*direction
    ds1 = ds(domain=mesh, subdomain_data=mf_1)
    l += inner(move,v)*ds1(2)
    deform = Function(W1)
    solve(a==l, deform , solver_parameters={"linear_solver":"mumps"})
    return deform

def lin_elasticity(mesh, mf_1, n, step, direction, cvt=None):
    # Rubber Beam Parameters
    rho = 950 # Density (kg/m^3)
    E = 0.1e9 # Youngs modulus (Pa)
    nu = 0.48 # Poisson ratio

    # Force due to gravity
    g = 9.81 # m/s^2
    f = Constant((0,0))

    # Elasticity parameters
    mu = E/(2.0*(1.0 + nu))
    lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

    # Stress computation
    def sigma(v):
        return 2.0*mu*sym(grad(v)) + lmbda*tr(sym(grad(v)))*Identity(len(v))

    # Create function space
    V = VectorFunctionSpace(mesh, "Lagrange", 1)


    # Stress on boundary
    x = SpatialCoordinate(mesh)
    dS_stress = Measure("ds", domain=mesh, subdomain_data=mf_1)
    if cvt is not None:
        move = step*n*direction+cvt
    else:
        move = step*n*direction

    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(sigma(u), grad(v))*dx
    L = inner(f, v)*dx # + inner(move, v)*dS_stress(2)

    bc = DirichletBC(V, move, mf_1, 2)

    # Create solution function
    u_fin = Function(V, name="deform")
    solve(a==L, u_fin, bcs=[bc], solver_parameters={"linear_solver": "mumps"})

    return u_fin


"""
Several helper functions for the linear and bilinear weak formulation
of the Stokes equation.
(see https://arxiv.org/pdf/1206.1933.pdf p. 5 eq 3.1-3.2)
"""
def b_h(v, q, n):
    return -div(v)*q*dX + jump(v, n)*avg(q)*dI

def l_h(v, q, f):
    return inner(f, v)*dX

def s_O(v, w):
    return inner(jump(grad(v)), jump(grad(w)))*dO

def s_C(v, q, w, r, h):
    return h*h*inner(-div(grad(v)) + grad(q), -div(grad(w)) - grad(r))*dC + \
        h("+")*h("+")*inner(-div(grad(v("+"))) + grad(q("+")),
                            -div(grad(w("+"))) - grad(r("+")))*dO
def l_C(v, q, f, h):
    return h*h*inner(f, -div(grad(v)) - grad(q))*dC\
        +  h("+")*h("+")*inner(f("+"), -div(grad(v("+"))) - grad(q("+")))*dO

def tensor_jump(v, n):
    return outer(v('+'), n('+')) + outer(v('-'), n('-'))

def a_h(v, w, n, h,  alpha=6.0): # FIXME: alpha=2,4. does not work,works: 6,10.
    return inner(grad(v), grad(w))*dX \
         - inner(avg(grad(v)), tensor_jump(w, n))*dI \
         - inner(avg(grad(w)), tensor_jump(v, n))*dI \
         + alpha/avg(h) * inner(jump(v), jump(w))*dI

"""
Returns the deformed mesh after getting gradient information
Input:
      omega - The inital mesh
      mf - MeshFunction with boundaries 
      gradient - gradient of functional in descent direction
      step - Step-lenght of scheme
      laplace - Type of deformation scheme
Output:
      omega - The deformed mesh
"""
def deform_mesh(multimesh_o, step, forget=False, vfac=Constant(1),
                bfac=Constant(1),cvt=None):
    if forget:
        multimesh = MultiMesh()
        multimesh.add(Mesh(multimesh_o.part(0)))
        multimesh.add(Mesh(multimesh_o.part(1)))
        multimesh.build()

    else:
        multimesh = multimesh_o

    V2 = VectorElement("CG", triangle, 2)
    Vdrag0 = FunctionSpace(multimesh.part(0), V2)
    Vdrag1 = FunctionSpace(multimesh.part(1), V2)
    Vdrag = MultiMeshFunctionSpace(multimesh, V2)

    u, p = StokesSolve(multimesh)
    gradient = functional_gradient(u, multimesh, Vol0, bx0, by0,vfac=vfac,
                                   bfac=bfac)
    direction = -gradient
    mf_0, mf_1 = load_facet_function(multimesh)
    
    from femorph import VolumeNormal
    normal = VolumeNormal(multimesh.part(1), [0], mf_1)

    deform_time = -time.time()
    # w1 =  Laplacian(multimesh.part(1), mf_1, normal, step, direction,
    #               alpha=5e-1)
    w1 = lin_elasticity(multimesh.part(1), mf_1, normal, step, direction
                        ,cvt=cvt)
    ALE.move(multimesh.part(1), w1)

    deform_time += time.time()
    print("Deformation time: %.2e" % deform_time)
    multimesh.build()
    return multimesh, w1

"""
Returns volume, baricenter of of a multimesh
Input:
      omega - The inital multimeshmesh
Output:
      Vol - Volume of multimesh
      bx,by - baricenter of multimesh
"""
def geometric_quantities(multimesh):
    mesh_time = -time.time()
    multimesh.build()
    multimesh.auto_cover(0, Point(0.5,0.5))

    mesh_time += time.time()
    print("Build and cover %.2e" % (mesh_time))
    x = SpatialCoordinate(multimesh)
    V = MultiMeshFunctionSpace(multimesh, "CG", 1)
    x0 = interpolate(Expression("x[0]", degree=1), V)
    x1 = interpolate(Expression("x[1]", degree=1), V)
    VolOmega = assemble_multimesh(1*dx(domain=multimesh)+1*dC(domain=multimesh))
    Vol = Constant(1 - VolOmega)
    # FIXME, Should be able to do something like
    # assemble_multimesh(x[0]*dx(domain=multimesh))
    bx = Constant((Constant(1./2)-assemble_multimesh(x0*dX))/Vol)
    by = Constant((Constant(1./2)-assemble_multimesh(x1*dX))/Vol)
    return Vol, bx, by

"""
Returns the functional value for a given mesh
Input:
      u - Solution of Stokes eq
      omega - MultiMesh for stokes-solution
      Vol0, bx0, by0 - Goal quantites for volume and baricenter
      vfac, bfac - Penalty parameters for volume and baricenter
Output:
     J - Functional value at given point in space
"""
def functional(u, multimesh, Vol0, bx0, by0, vfac=Constant(1e5), bfac=Constant(1e3)):
    Vol, bx, by = geometric_quantities(multimesh)
    vol_off = Vol-Vol0
    bx_off, by_off = bx-bx0, by-by0
    J_s = assemble_multimesh(inner(grad(u),grad(u))*dX)
    J_v = vfac*vol_off*vol_off
    J_c = bfac*(bx_off*bx_off+by_off*by_off)
    J = J_s + J_v + J_c
    return J

def functional_gradient(u, multimesh, Vol0, bx0, by0, vfac=Constant(1e5),
                        bfac=Constant(1e3)):
    x = SpatialCoordinate(multimesh.part(1))
    u_1 = u.part(1, deepcopy=True)
    Vol, bx, by = geometric_quantities(multimesh)
    vol_off = Vol-Vol0
    bx_off, by_off = bx-bx0, by-by0
    stokes = -inner(grad(u_1), grad(u_1))
    vol = - Constant(2.)*vfac*vol_off
    bar = Constant(2.)*bfac/Vol*((bx-x[0])*bx_off + (by-x[1])*by_off)
    return stokes+vol+bar

def splitMMF(w,multimesh):
    V2 = VectorElement("CG", triangle, 2)
    S1 = FiniteElement("CG", triangle, 1)
    V0 = FunctionSpace(multimesh.part(0), V2)
    V1 = FunctionSpace(multimesh.part(1), V2)
    Vmm = MultiMeshFunctionSpace(multimesh, V2)
    P0 = FunctionSpace(multimesh.part(0), S1)
    P1 = FunctionSpace(multimesh.part(1), S1)
    Pmm = MultiMeshFunctionSpace(multimesh, S1)
    u0, p0 = w.part(0, deepcopy=True).split()
    u1, p1 = w.part(1, deepcopy=True).split()
    umm = MultiMeshFunction(Vmm)
    pmm = MultiMeshFunction(Pmm)
    umm.assign_part(0, interpolate(u0,V0))
    umm.assign_part(1, interpolate(u1,V1))
    pmm.assign_part(0, interpolate(p0,P0))
    pmm.assign_part(1, interpolate(p1,P1))
    return umm,pmm

if __name__ == "__main__":
    import sys
    import time
    start = time.time()
    
    # Get mesh and meshfunction from file
    meshes =[]
    multimesh = MultiMesh()
    for i in range(2):
        mesh_i = Mesh()
        with XDMFFile("meshes/multimesh_%d.xdmf" %i) as infile:
            infile.read(mesh_i)
        meshes.append(mesh_i)
        multimesh.add(mesh_i)
    multimesh.build()
            

    V2 = VectorElement("CG", triangle, 2)
    S1 = FiniteElement("CG", triangle, 1)
    TH = V2 * S1
    VQ = MultiMeshFunctionSpace(multimesh, TH)
    us,ps = [],[]
    for i in range(multimesh.num_parts()):
        us.append(XDMFFile("output/u%d.xdmf" %i))
        ps.append(XDMFFile("output/p%d.xdmf" %i))
    out = {"VQ": VQ,"w":MultiMeshFunction(VQ)}

    # Compute original volume and baricenter of obstacle
    Vol0, bx0, by0 = geometric_quantities(multimesh)

    it, max_it, mq_tol, mq, func_old = 0, 50, 0.05, 0.5, 1e2
    v_o = Constant(0)
    vfac, bfac = Constant(1e5), Constant(1e3)
    red_tol = 1e-5 # 8e-5 standard red_tol
    start_stp = 1e-4
    stp_min,stp_max = 5e-9, 1e-1

    sub_problem_it, Js, dJs, Vol_off, Bx_off, By_off, MQ = ([] for _ in range(7))
    while (it<max_it):
        time_it = time.time()
        if float(v_o)!=float(vfac):
            if it != 0:
                stp_min/=2.
                stp_max*=2.
            sub_problem_it.append(it)
            print("Linesearch initialized")
            print(start_stp, stp_max)
            search = moola.linesearch.ArmijoLineSearch(start_stp=start_stp,
                                                       stpmax=stp_max,
                                                       stpmin=stp_min)
            v_o = vfac
        u, p = StokesSolve(multimesh, out=out)
        for dof in u.function_space().dofmap().inactive_dofs(multimesh,0):
            u.vector()[dof]=np.nan
        for dof in p.function_space().dofmap().inactive_dofs(multimesh,0):
            p.vector()[dof]=np.nan
        for i in range(multimesh.num_parts()):
            us[i].write(u.part(i, deepcopy=True), float(it))
            ps[i].write(p.part(i, deepcopy=True), float(it))
        J = functional(u, multimesh, Vol0, bx0, by0, vfac=vfac, bfac=bfac)
        gradient = functional_gradient(u, multimesh, Vol0, bx0, by0,
                                       vfac=vfac, bfac=bfac)
        Js.append(float(J))
        
        def phi(step):
            step = Constant(step)
            print("Armijo-functional step %.3e" % step)
            multimesh_s, w1_s = deform_mesh(multimesh, step, forget=True,
                                            vfac=vfac, bfac=bfac)
            u_s, p_s = StokesSolve(multimesh_s)
            J = functional(u_s, multimesh_s, Vol0, bx0, by0, vfac=vfac,
                           bfac=bfac)
            return float(J)
        def phi_dphi0():
            [mf_0, mf_1] = load_facet_function(multimesh)
            ds1 = ds(domain=multimesh.part(1), subdomain_data=mf_1)
            from femorph import VolumeNormal
            normal = VolumeNormal(multimesh.part(1), [0], mf_1)
            dJdstep = assemble(inner(gradient*normal, -gradient*normal)*ds1(2))
            return float(J), dJdstep
        def taylor():
            steps = [10**-(i+5) for i in range(7)]
            for step in steps:
                J_new = phi(step)
                print(float((J_new-J)/(step)), "FD")
            print(phi_dphi0()[1], "gradient")   
        # taylor()

        # try:
        print('*' *14 + "Armijo linesearch" + "*"*14)
        time_search = -time.time()
        step = Constant(search.search(phi, None, phi_dphi0()))
        time_search += time.time()
        print("Linesearch time: %.2e" %(time_search))

        # except:
        #     # sub_problem_it.append(it)
        #     print("Reached minimal step-length, increasing penalty")

        #     vfac = Constant(2*vfac)
        #     bfac = Constant(2*bfac)
        #     it+=1
        #     continue
        tmp_diff = func_old-J
        func_old = J
        Vol, bx,by = geometric_quantities(multimesh)
        vol_off, bx_off, by_off = Vol-Vol0, bx-bx0, by-by0
        tmp_grad_norm = phi_dphi0()[1]
        dJs.append(np.sqrt(np.abs(tmp_grad_norm)))
        Vol_off.append(float(vol_off))
        Bx_off.append(float(bx_off))
        By_off.append(float(by_off))
        mq_pre = min(MeshQuality.radius_ratio_min_max(multimesh.part(0))[0],
                 MeshQuality.radius_ratio_min_max(multimesh.part(1))[0])

        print('*'*5 + "CVT-smoothing" + '*'*5)
        from femorph.Legacy import DiscreteMeshRepair
        from femorph import VolumeNormal
        mf_0,mf_1 = load_facet_function(multimesh)
        cvt_time = -time.time()
        n_ = VolumeNormal(multimesh.part(1),[0], mf_1)
        (L1B, L2B,fix_deform) = DiscreteMeshRepair(multimesh.part(1), mf_1,
                                                   SmoothVolume=False,
                                                   SmoothBoundary=True,
                                                   Tangential=True,
                                                   VertexNormal=n_,
                                                   MaxIter = 10, Step = 1.0,
                                                   Stop = 1e-6,
                                                   FixPlanes = [[1,0.5,1e-3]],Vis=False)
        cvt_time += time.time()
        print("CVT-time: %.2e" % cvt_time) 
        print("Updating domain")
        
        
        print("Updating domain")
        multimesh, w1 = deform_mesh(multimesh, step, vfac=vfac, bfac=bfac,
                                    cvt=fix_deform)
        mq = min(MeshQuality.radius_ratio_min_max(multimesh.part(0))[0],
                 MeshQuality.radius_ratio_min_max(multimesh.part(1))[0])
        MQ.append(mq)
        print("*"*45)
        print("Iteration: %d" %it)
        print("J:       %.5f" % (J))
        print("Relative decrease : %.3e" % (-tmp_diff/J))
        print("V_off, B_off:  %.2f%% %.2f%% %.2f%%"\
            %(float((vol_off)/(Vol0)*100),float(bx_off/bx0*100),
              float(by_off/by0*100)))
        print("Gradient norm: %.5e" %np.sqrt(np.abs(tmp_grad_norm)))
        print("Mesh quality %.2f" % mq_pre)

        if abs(float(tmp_diff/J)) < red_tol:
            vfac = Constant(2*vfac)
            bfac = Constant(2*bfac)
            print("Reduction tolerance reached (%.2e), increasing penalty" %red_tol)
            if float(vfac) > 2**4*1e5:
                print("Reached max penalty parameter, Volume offset: %.2e%%" %(float((vol_off)/(Vol0)*100)))
                break
        it += 1
        print("Iteration-time: %d" %(time.time()-time_it)) 
        if mq < 0.05:
            plot(multimesh.part(1))
            plt.show()
            raise ValueError("Mesh Degenerated")

    np.savez("output/Optimization_results.npz", J=Js,dJ=dJs, Vol=Vol_off, Bx=Bx_off, By=By_off, MQ=MQ, subproblem=sub_problem_it)
    print("Endtime %3d" %(time.time()-start))
