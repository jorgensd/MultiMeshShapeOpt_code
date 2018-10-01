
import numpy
import dolfin
from dolfin import *
import pyipopt
from IPython import embed
import sympy
from matplotlib import pylab as plt
# from Cables_3 import eval_J, eval_dJ, update_mesh, T, cable_positions, num_cables, cable_scales
# from Cables_5 import eval_J, eval_dJ, update_mesh, T, cable_positions, num_cables, cable_scales
from Cable_isosceles import  eval_J, eval_dJ, update_mesh, T, cable_positions, num_cables, cable_scales 
set_log_level(LogLevel.ERROR)

# Initial cable information
scale = 100
outer_cable_mesh_radius = scale*0.012  
inner_cable_mesh_radius = scale*0.003*cable_scales  # This includes mesh halo
distance_from_outer_cable_mesh = scale*0.001
max_radius =  outer_cable_mesh_radius - inner_cable_mesh_radius  \
              - distance_from_outer_cable_mesh


# Constraint informatuio
g_scale = 1e+5  # Scaling coefficient for inequality constraint
def sympy_g(num_cables):
    """ Create constraint for N cables """
    x = list(sympy.symbols("x0:%d" %num_cables))
    y = list(sympy.symbols("y0:%d" %num_cables))
    g,z= [],[]
    for i in range(num_cables):
        g.append(x[i]**2 + y[i]**2 - max_radius[i]**2)
    for i in range(num_cables):
        z.append(x[i])
        z.append(y[i])
        for j in range(i+1,num_cables):
            xi, yi = x[i],y[i]
            xj, yj = x[j],y[j]
            int_radius = inner_cable_mesh_radius[i]+inner_cable_mesh_radius[j]\
                         +distance_from_outer_cable_mesh
            # FIXME: Max_int_radius should change with each cable radius
            g.append(int_radius**2-(xi-xj)**2-(yi-yj)**2)
    g = sympy.Matrix(g)
    z = sympy.Matrix(z)
    return g,z

# Creating constraints with sympy
g,z = sympy_g(num_cables)

def replace_sympy_g(positions):
    """ Create numpy array of constraint at current positions """
    g_numpy = g.subs([(z[i], positions[i]) for i in range(2*num_cables)])
    return g_scale*numpy.array(g_numpy).astype(float)

def sympy_jac_g():
    """ Creates jacobian for constraints """
    jac_g = g.jacobian(z)
    return jac_g

jac_g = sympy_jac_g()

def replace_sympy_jac_g(positions):
    """ Create numpy array of jacobian of constraint at current positions """    
    numpy_jac_g = jac_g.subs([(z[i], positions[i]) for i in range(2*num_cables)])
    return g_scale*numpy.array(numpy_jac_g).astype(float)



def eval_g(cable_positions):
    """ Inequality constraint, g(x) <= 0 """
    return replace_sympy_g(cable_positions)

def eval_jac_g(cable_positions, flag):
    """ The constraint Jacobian:
        flag = True  means 'tell me the sparsity pattern';
        flag = False means 'give me the Jacobian'."""
    if flag:
        # pass in a dense matrix (could be optimised, but we don't care for this small problem).
        nvar = len(cable_positions)
        n_cables = int(nvar/2)
        ncon = int(n_cables + n_cables*(n_cables-1)/2)
        rows = []
        for i in range(ncon):
            rows += [i] * nvar
        cols = list(range(nvar)) * ncon
        return (numpy.array(rows), numpy.array(cols))
    else:
        return replace_sympy_jac_g(cable_positions)

def main():
    cable_positions_numpy = numpy.array(cable_positions)
    
    print( cable_positions)
    def eval_robust_J(cable_positions):
        if max(eval_g(cable_positions)) > 0:
            print( "Warning: Inequality constraint is not satisfied")
            return numpy.inf
        # print eval_g(cable_positions).T
        try:
            j = eval_J(cable_positions)
            # c1 = numpy.array([cable_positions[0],cable_positions[1]])
            # c2 = numpy.array([cable_positions[2],cable_positions[3]])
            # c3 = numpy.array([cable_positions[4],cable_positions[5]])
            # print("Angles between the three cables: %.2f, %.2f, %.2f"
            #       % (numpy.arccos(numpy.dot(c1-c2,c1-c3)/
            #                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
            #                                   *numpy.dot(c1-c3,c1-c3))))/(2*numpy.pi)*360,
            #          numpy.arccos(numpy.dot(c1-c2,c3-c2)/
            #                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
            #                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360,
            #          numpy.arccos(numpy.dot(c3-c1,c3-c2)/
            #                       (numpy.sqrt(numpy.dot(c3-c1,c3-c1)
            #                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360))
        except:
            print("ERROR: Forward model failed, returning infinity")
            j = numpy.inf
            exit(1)

        #print( "J = ", j)
        return j

    def eval_robust_dJ(cable_positions):
        # if max(eval_g(cable_positions)) > 0:
        #     print( "Warning: Inequality constraint is not satisfied")
            # return numpy.inf
        dj = eval_dJ(cable_positions)
        print("|dj| = ", numpy.dot(dj, dj)**0.5)
        return dj

    def test_J(cable_positions):
        """ Test of optimization without state eq"""
        j = 0.5*numpy.dot(cable_positions, cable_positions)
        if max(eval_g(cable_positions)) > 0:
            print("Warning: Inequality constraint is not satisfied: ")
            print(eval_g(cable_positions))
            print(cable_positions)
            print("*"*25)
        print("J = ", j)
        return j

    def test_dJ(cable_positions):
        """ Test of optimization without state eq"""
        dj = numpy.array(cable_positions)
        print("|dj| = ", numpy.dot(dj, dj)**0.5)
        return dj

    # Check that input data match input variables
    nvar = int(2*num_cables) # Number of controls
    ncon = int(num_cables + num_cables*(num_cables-1)/2) # Number of inequality constraints 

    # Create the NLP model
    #pyipopt.set_loglevel(2)         # Set verbosity
    nlp = pyipopt.create(
        nvar,                                          # Number of controls
        -numpy.inf*numpy.ones(nvar,dtype=float),   # Lower bounds of controls
        numpy.inf*numpy.ones(nvar, dtype=float),   # Upper bounds of controls
        ncon,                                          # Number of inequality constraints
        -numpy.inf*numpy.ones(ncon, dtype=float),   # Lower bounds of inequality constraints
        numpy.zeros(ncon, dtype=float),            # Upper bounds of inequality constraints
        nvar*ncon,                                     # Number of nonzeros in the constraint Jacobian
        0,                                             # Number of nonzeros in the Hessian
        # lambda pos: eval_robust_J(pos),                # Objective evaluation
        # lambda pos: eval_robust_dJ(pos),               # Objective gradient evaluation
        lambda pos: -test_J(pos),                      # Objective evaluation
        lambda pos: -test_dJ(pos),                     # Objective gradient evaluation
        eval_g,                                        # Constraint evaluation
        eval_jac_g,                                    # Constraint Jacobian evaluation
    )
    nlp.num_option('obj_scaling_factor',1e-2)#1e-7)
    nlp.int_option('max_iter', 100)
    nlp.int_option('print_level', 5)
    nlp.num_option('tol', 1e-4)
    nlp.str_option('mu_strategy',"adaptive")
    nlp.num_option("mu_max", 1)
    nlp.num_option('mu_init', 1)#1e25)
    nlp.num_option('bound_relax_factor', 0) # So it does not violate the boundary constraint
    nlp.num_option('acceptable_tol', 1e-5) 
    # Solve the optimisation problem
    opt_cable_pos = nlp.solve(cable_positions_numpy)[0]
    nlp.close()
    
    # Report the results
    d = opt_cable_pos.reshape(-1, 2)
    for i in range(num_cables):
        max_radius = outer_cable_mesh_radius - inner_cable_mesh_radius - distance_from_outer_cable_mesh
        print("Cable %d distance from center: %.3e" %(i,numpy.sqrt(d[i,0]**2+d[i,1]**2)))
        for j in range(i+1, num_cables):
            print("Distance between Cable %d and %d: %.3e" %(i,j, numpy.sqrt((d[i,0]-d[j,0])**2)
                                                                             +(d[i,1]-d[j,1])**2))
    print("Cable_positions:")
    for i in range(num_cables):
        print("%10.3e, %10.3e" %(d[i,0],d[i,1]))
    update_mesh(opt_cable_pos)

    # dolfin.plot(update_mesh(opt_cable_pos))
    TempFiles = [dolfin.File("output/T_ipopt_part%d.pvd" %i) for i in range(num_cables+1)]
    eval_J(cable_positions)
    for i in range(num_cables+1):
        TempFiles[i] << T.part(i)
    
    eval_J(opt_cable_pos)    
    for i in range(num_cables+1):
        TempFiles[i] << T.part(i)
    c1 = d[0,:]
    c2 = d[1,:]
    c3 = d[2,:]
    print(numpy.arccos(numpy.dot(c1-c2,c1-c3)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c1-c3,c1-c3))))/(2*numpy.pi)*360)
    print(numpy.arccos(numpy.dot(c1-c2,c3-c2)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)
    print(numpy.arccos(numpy.dot(c3-c1,c3-c2)/
                       (numpy.sqrt(numpy.dot(c3-c1,c3-c1)
                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)



if __name__ == '__main__':
    main()
