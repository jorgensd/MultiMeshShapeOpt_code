from dolfin import *
import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
# Physical parameters
num_cables = 3       # Number of inner cables
lmb_metal = 205.      # Heat coefficient inside the metal (aluminium)
lmb_insulation = 0.03  # Heat coefficient of plastic insulating material
lmb_air =0.15  # Heat coefficient of brick (filling material in cable)

# Load cables
cable_meshes = []
cable_subdomains = []
cable_facets = []

def projection(cable_positions, max_radius):
    """
    Projecting the inner cables into the larger cable if they are moved outside
    of the valid domain, returns the new cable positions if projected, as well
    as number of violations
    """
    cable_positions_tuples = cable_positions.reshape(-1, 2)
    icount = 0
    for i, (x, y) in enumerate(cable_positions_tuples):
        # A projected gradient method, see Algorithm 2.3, The
        # optimality condition in Corollary 1.2, Hinze and Ulbrich 2009
        if x**2 + y**2 > max_radius[i]**2:
            icount+=1
            alpha = max_radius[i]/sqrt(x**2+y**2)
            cable_positions[2*i] = alpha*cable_positions[2*i]
            cable_positions[2*i+1] = alpha*cable_positions[2*i+1]
            assert (cable_positions[2*i]**2 + cable_positions[2*i+1]**2
                    <= max_radius[i]**2 + 3*DOLFIN_EPS)
    if icount>0:
        print( "%d cables outside domain, projecting into domain " %icount)
    return cable_positions, icount


# Initial cable positions
# Correct angles
# rad = 0.36 # in [0.35-1]
# c1 = rad*numpy.array([numpy.cos(0), numpy.sin(0)])
# c2 = rad*numpy.array([numpy.cos(-2*numpy.pi/3),numpy.sin(-2*numpy.pi/3)])
# c3 = rad*numpy.array([numpy.cos(2*numpy.pi/3), numpy.sin(2*numpy.pi/3)])

# Wrong angles close to center
c1 = numpy.array([-0.2, 0.4])
c2 = numpy.array([-0.35, -0.3])
c3 = numpy.array([0.3,-0.4])



cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
scale = 100
outer_cable_mesh_radius = scale*0.012  
inner_cable_mesh_radius = scale*0.003*numpy.array([1,1,1])
distance_from_outer_cable_mesh = scale*0.0001
max_radius = (outer_cable_mesh_radius - inner_cable_mesh_radius
            - distance_from_outer_cable_mesh)
#cable_positions, d = projection(cable_positions, max_radius)
cable_scales = numpy.array([1, 1,1])

def compute_angles(cable_positions):
    c1 = cable_positions[0:2]
    c2 = cable_positions[2:4]
    c3 = cable_positions[4:6]
    a1,a2,a3 = (numpy.arccos(numpy.dot(c1-c2,c1-c3)/
                             (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                         *numpy.dot(c1-c3,c1-c3))))/(2*numpy.pi)*360,
                numpy.arccos(numpy.dot(c1-c2,c3-c2)/
                             (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                         *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360,
                numpy.arccos(numpy.dot(c3-c1,c3-c2)/
                             (numpy.sqrt(numpy.dot(c3-c1,c3-c1)
                                        *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)

    print("Angles between the three cables: %.2f, %.2f, %.2f" %(a1,a2,a3))
compute_angles(cable_positions)

# Background mesh
mesh_0 = Mesh("meshes/cable.xml")
facets_0 = MeshFunction("size_t", mesh_0, "meshes/cable_facet_region.xml")

# Move origo-centered meshes to initial posiion
for i in range(num_cables):
    cable = Mesh("meshes/inner_cable_halo.xml")
    cable.coordinates()[:,:]*=cable_scales[i]
    cable.translate(Point(cable_positions[2*i],cable_positions[2*i+1]))
    cable_meshes.append(cable)

    subdomains = MeshFunction("size_t", cable, "meshes/inner_cable_halo_physical_region.xml")
    cable_subdomains.append(subdomains)

    facets = MeshFunction("size_t", cable, "meshes/inner_cable_halo_facet_region.xml")
    cable_facets.append(facets)

# Create multimesh
multimesh = MultiMesh()
multimesh.add(mesh_0)
for cable in cable_meshes:
    multimesh.add(cable)
multimesh.build()

# Create function space for temperature
V = MultiMeshFunctionSpace(multimesh, "CG", 1)
W = MultiMeshFunctionSpace(multimesh, "DG", 0)
V0 = FunctionSpace(mesh_0, "CG", 1)
W0 = FunctionSpace(mesh_0, "DG", 0)
Vcables = [FunctionSpace(cable, "CG", 1) for cable in cable_meshes]

# Construct forcing function
# Note: we assume that all cable meshes have the same data structure
f = MultiMeshFunction(W)
f.assign_part(0, project(Constant(0.0), W0))
X = FunctionSpace(cable_meshes[0], "DG", 0)
fx = Function(X)
# cnorm = mpl.colors.Normalize(vmin=0, vmax=100)

for i in range(num_cables):
    fx.vector()[:] = (14.5 < cable_subdomains[0].array() )*(5 + 10*(i>1))
    f.assign_part(i+1, fx)
    # plot(f.part(i+1, deepcopy=True), norm=cnorm)
# Construct heat coefficient function
# Note: we assume that all cable meshes have the same data structure
lmb = MultiMeshFunction(W)
lmb.assign_part(0, project(Constant(lmb_air), W0))
lmbx = Function(X)
lmbx.vector()[:] = ((cable_subdomains[0].array() == 13)*lmb_air + 
                    (cable_subdomains[0].array() == 14)*lmb_insulation +
                    (cable_subdomains[0].array() == 15)*lmb_metal)

for i in range(num_cables):
    lmb.assign_part(i+1, lmbx)
# Define trial and test functions and right-hand side
T = MultiMeshFunction(V, name="Temperature")
adjT = MultiMeshFunction(V, name="Adjoint")

# Define facet normal and mesh size
n = FacetNormal(multimesh)
h = 2.0*Circumradius(multimesh)
h = (h('+') + h('-')) / 2

# Define multimesh parameters
alpha = 4.0
beta = 4.0

# Define objective
q = 3.0#5.0
obj = (1.0/q)*(pow(abs(T), q))*dX
objdT = T*pow(abs(T), q-2)

# Define coefficients for the state equation
T_amb = Constant(3.2)                # ambient Temperature
c = Constant(0.04)#Constant(0.1)
def alpha_heat_transfer(T):          # heat transfer coefficient at ex bc
    return Constant(1.0)

""" The Riez representer of the shape-surface gradient """
def WeakCableShapeGradSurf(T, adjT, lmb, c, f, n):
    # ("-") are exterior quantitites, [ ]_{+-}=-jump( )
    grad_tau_adjT = grad(adjT("-"))-dot(n,grad(adjT("-")))*n 
    grad_tau_T = grad(T("-"))-dot(n,grad(T("-")))*n
    dJ = inner(grad_tau_T, grad_tau_adjT)*jump(lmb)\
         - adjT("-")*(jump(f)# +jump(c)*adjT("-")*T("-")
         )\
         - lmb("-")*dot(n,grad(adjT("-")))*jump(dot(n, grad(T)))
    # dJ = jump(inner(grad(adjT)-dot(n,grad(adjT))*n,
    #                 grad(T)-dot(n,grad(T))*n)*lmb) \
    #                 -jump(adjT*f)-jump(c*T*adjT) \
    #                 -jump(lmb*dot(n,grad(adjT))*dot(n, grad(T)))
    return -dJ


""" Translate all new_cables to a new center """
def update_mesh(cable_positions,angles=False):
    cable_positions = cable_positions.reshape(-1, 2)
    for i, ((newx, newy), cable_mesh, cable_facet) in enumerate(
            zip(cable_positions, cable_meshes, cable_facets)):
        # Move mesh to new positions
        dSc = Measure("dS", subdomain_data=cable_facet, subdomain_id=17)
        area = assemble(Constant(1)*dx(domain=cable_mesh))

        oldx = assemble(Expression("x[0]", degree=1)*dx(domain=cable_mesh))/area
        oldy = assemble(Expression("x[1]", degree=1)*dx(domain=cable_mesh))/area
        diff = [-oldx+newx, -oldy+newy]
        # Only rebuild the mesh if there is actually any movement of the cables
        cable_mesh.translate(Point(-oldx+newx, -oldy+newy))
        # Finding angle between radial vector and perturbation direction
        pertb = numpy.array([-oldx+newx,-oldy+newy])\
                /numpy.sqrt((-oldx+newx)**2+(-oldy+newy)**2)
        i_r = numpy.array([oldx, oldy])/numpy.sqrt(oldx**2+oldy**2)
    #     if angles:
    #          print( "Angles Cable %d, i_r: %.2f"
    #                 %(i, numpy.arccos(numpy.dot(i_r, pertb)
    #                                   /numpy.sqrt(numpy.dot(i_r,i_r)
    #                                               *numpy.dot(pertb,pertb)))
    #                   *180/numpy.pi))
    # print("")
    multimesh.build()  # Rebuild the multimesh
    return multimesh

""" Evaluate the functional with given cable_positions"""
def eval_J(cable_positions, angles=False):
    # Update mesh
    update_mesh(cable_positions, angles)

    v = TestFunction(V)
    Ttmp = TrialFunction(V)
    constraint = inner(lmb*grad(Ttmp), grad(v))*dX \
             -f*v*dX -c*v*Ttmp*dX
    constraint += alpha_heat_transfer(Ttmp)*(Ttmp-T_amb)*v*ds
    constraint += - inner(avg(lmb*grad(Ttmp)), jump(v, n))*dI \
                  - inner(avg(lmb*grad(v)), jump(Ttmp, n))*dI \
                  + alpha/h*jump(Ttmp)*jump(v)*dI   \
                  + beta*lmb*inner(jump(grad(Ttmp)), jump(grad(v)))*dO
    A, b= assemble_multimesh(lhs(constraint)), assemble_multimesh(rhs(constraint))
    T.vector()[:]=0
    V.lock_inactive_dofs(A, b)
    solve(A, T.vector(), b, 'lu')
    
    return assemble_multimesh(obj)

# Evaluate the shape gradient
def eval_dJ(cable_positions):
    # Update mesh
    eval_J(cable_positions)
    dJ = []
    # Solve adjoint equation
    adj = TrialFunction(V)
    v = TestFunction(V)
    constraint = inner(lmb*grad(adj), grad(v))*dX -c*v*adj*dX
    # FIXME: Add derivative of alpha in ext bc constraint,only works for alpha=1
    constraint += adj*1*v*ds# alpha_heat_transfer(T)*v*ds
    constraint += - inner(avg(lmb*grad(adj)), jump(v, n))*dI \
                  - inner(avg(lmb*grad(v)), jump(adj, n))*dI \
                  + alpha/h*jump(adj)*jump(v)*dI   \
                  + beta*lmb*inner(jump(grad(adj)), jump(grad(v)))*dO
    constraint += objdT*v*dX
    A, b= assemble_multimesh(lhs(constraint)), assemble_multimesh(rhs(constraint))
    V.lock_inactive_dofs(A, b)
    solve(A, adjT.vector(), b, 'lu')

    for i, (cable_mesh, cable_facet,cable_subdomain) in enumerate(
            zip(cable_meshes, cable_facets,cable_subdomains)):
        T_cable = T.part(i+1, deepcopy=True)
        adjT_cable = adjT.part(i+1, deepcopy=True)
        lmb_cable = lmb.part(i+1, deepcopy=True)
        f_cable = f.part(i+1, deepcopy=True)

        from femorph import VolumeNormal
        normal = VolumeNormal(cable_mesh, [0], cable_facet, [16,17])
        dJ_Surf = WeakCableShapeGradSurf(T_cable, adjT_cable,
                                         lmb_cable, c, f_cable, n=normal)
        dSc1 = Measure("dS", subdomain_data=cable_facet, subdomain_id=16)
        dSc2 = Measure("dS", subdomain_data=cable_facet, subdomain_id=17)
        gradx = assemble(normal[0]*dJ_Surf*dSc1
                         +normal[0]*dJ_Surf*dSc2
                         + Constant(0)*dx(domain=cable_mesh,
                                         subdomain_data=cable_subdomain))
        grady = assemble(normal[1]*dJ_Surf*dSc1
                         +normal[1]*dJ_Surf*dSc2
                         + Constant(0)*dx(domain=cable_mesh,
                                         subdomain_data=cable_subdomain))
       
        dJ.append(gradx)
        dJ.append(grady)
    return numpy.array(dJ)

def main():
    # Steepest decent is not reliable
    TFiles = [File("output/isoT_part_%d.pvd" %(i)) for i in range(num_cables+1)]
    # FFiles = [File("output/isof-part_%d.pvd" %(i)) for i in range(num_cables+1)]
    # LFiles = [File("output/isolmb-part_%d.pvd" %(i)) for i in range(num_cables+1)]

    # Implementation of a steepest descent method with projection
    # for the inequality constraints
    scale = 100
    outer_cable_mesh_radius = scale*0.012
    r_met, r_iso = 0.2*cable_scales, 0.255*cable_scales
    inner_cable_mesh_radius = scale*0.003*cable_scales  # This includes mesh hal
    distance_from_outer_cable_mesh = scale*0.0001
    max_radius = (outer_cable_mesh_radius - inner_cable_mesh_radius
                  - distance_from_outer_cable_mesh)

    #Optimization loop
    old_cable_positions = numpy.inf*numpy.ones(cable_positions.shape)
    local_cable_positions = cable_positions.copy()
    old_loc = local_cable_positions.copy()
    it = 0
    max_iter = 250
    tol = 1e-8
    js = [eval_J(local_cable_positions, True)]
    step = 1
    min_stp = 1e-6
    max_stp = 1
    last_step = step
    prog = 1e10
    rel_tol = 1e-6
    print("Starting optimization")
    #from plot_3_cables import plot_T #,plot_lambda
    T_max = numpy.max(T.vector().get_local()) # Max temperature at initial configuration
    while it < max_iter:
        if prog ==0:
            break
        # Save for each iteration
        # k-th cable position
        old_cable_positions = local_cable_positions.copy()
        gradient = eval_dJ(old_cable_positions)
        if it >0:
            js.append(Jtmp)
            j = Jtmp
        else:
            j = js[0]
        print( "-"*30, "Iter {}".format(it), "-"*30)
        print( "Functional value: {}".format(j))
        compute_angles(old_cable_positions)
        for i in range(num_cables+1):
            TFiles[i] << T.part(i)
            # FFiles[i] << f.part(i)
            # LFiles[i] << lmb.part(i)
        print("Relative decrease %.2e" %( numpy.abs(js[it]-js[it-1])/js[it]))
        if it == 0:
            #plot_T(T, T_max,old_cable_positions, r_met, r_iso, "Initial_temperature_iso")
            pass
        elif numpy.abs(js[it]-js[it-1])/js[it] < rel_tol:
            print("Decrease less than reduction tolerance, optimized solution found")
            break
            

        # Projected Armijo rule
        # (section 2.2.2.1 in Hinze and Ulbrich, 2009)
        proj_armijo=False
        grad_de = 1e-4 # default f-tol from moola
        while not proj_armijo and step>min_stp:
            # k+1-th cable position
            tmp_old = old_cable_positions.copy()
            assert((old_loc==old_cable_positions).all())
            local_cable_positions = old_cable_positions - step*gradient
            # Projected k+1-th cable position
            #print("Angle before projection")
            print(abs(step*gradient), "Relative movement")
            update_mesh(local_cable_positions, True)
            local_cable_positions, violate = projection(local_cable_positions,
                                                        max_radius)
            # k+1-th functional value
            Jtmp = eval_J(local_cable_positions, False)
            # Evaluation of projected Armijo rule
            if not (Jtmp-j< -grad_de/step*numpy.linalg.norm(
                    local_cable_positions-old_cable_positions)**2):
                print( "Projected Armijo not fulfilled, halving step-length",
                       " to %.2e" %(0.5*step), Jtmp-j, ">=", -grad_de/step*numpy.linalg.norm(local_cable_positions-old_cable_positions)**2)
                prog = numpy.linalg.norm(old_cable_positions
                                         - local_cable_positions)
                print( "Armijo Progress ",prog)
                step = step*0.5
            else:
                # Next iteration if satisfied
                proj_armijo=True
                print( "Projected armijo fulfilled with step %.3e" % step)
                prog = numpy.linalg.norm(old_cable_positions
                                         - local_cable_positions)
                print( "Progress ",prog)
            update_mesh(old_cable_positions)
            assert((tmp_old == old_cable_positions).all())

        # End if reached minimal step-length
        if step <= min_stp:
            print( "Minimum step-length occurred, found minima")
            break
        elif step>=max_stp:
            step = max_stp
            last_step=step
        else:
            step*=2
            last_step=step
        old_loc = local_cable_positions.copy()
        it += 1
    # Save final temperature
    plot_T(T, T_max,old_cable_positions, r_met, r_iso,"Final_temperature_iso")
    
    print( "Initial J: %.5e" % js[0])
    print( "Optimal J: %.5e" % js[-1])
    d = old_cable_positions.reshape(-1, 2)
    for i in range(num_cables):
        max_radius = outer_cable_mesh_radius - inner_cable_mesh_radius \
                     - distance_from_outer_cable_mesh
        print( "Cable %d distance from center: %.3e" %(i,numpy.sqrt(d[i,0]**2+d[i,1]**2)))
        for j in range(i+1, num_cables):
            print( "Distance between Cable %d and %d: %.3e" %(i,j, numpy.sqrt((d[i,0]-d[j,0])**2
                                                                             +(d[i,1]-d[j,1])**2)))
    print( "Cable_positions:")
    for i in range(num_cables):
        print( "%10.3e, %10.3e" %(d[i,0],d[i,1]))
    compute_angles(old_cable_positions)

if __name__ == "__main__":
    main()
 
