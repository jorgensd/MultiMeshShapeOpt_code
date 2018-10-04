from dolfin import *
from IPython import embed
import numpy
import pyipopt
import sympy


class IpoptAngle():
    
    def __init__(self, num_angles, J, dJ):
        self.num_angles = num_angles
        self.solve_init(J, dJ)

    def func_g(self, x, user_data=None):
        return empty

    def jac_g(self, x, flag, user_data=None):
        if flag:
            rows = numpy.array([], dtype=int)
            cols = numpy.array([], dtype=int)
            return (rows, cols)
        else:
            empty = numpy.array([], dtype=float)
            return empty

        
    def solve_init(self, eval_J, eval_dJ):
        nvar = self.num_angles
        low_var = -numpy.infty*numpy.ones(nvar,dtype=float)
        up_var = numpy.infty*numpy.ones(nvar, dtype=float)


        self.nlp = pyipopt.create(nvar,     # Number of controls
                                  low_var,  # Lower bounds for Control
                                  up_var,   # Upper bounds for Control
                                  0,     # Number of constraints
                                  numpy.array([], dtype=float),
                                  numpy.array([], dtype=float),
                                  0,        # Number of nonzeros in cons. Jac
                                  0,        # Number of nonzeros in cons. Hes
                                  lambda angle: eval_J(angle),  # Objective eval
                                  lambda angle: eval_dJ(angle), # Obj. grad eval
                                  self.func_g,
                                  self.jac_g)

    def solve(self, angle):
        return self.nlp.solve(angle)[0]


if __name__ == "__main__":
    from StokesSolver import *
    import numpy
    Points = [0.25,0.25,0.5,0.25,0.75,0.25,
              0.25,0.5, 0.5,0.5, 0.75,0.5,
              0.25,0.75,0.5,0.75,0.75,0.75]
    points = [Point(Points[2*i],Points[2*i+1])
              for i in range(int(len(Points)/2))]
    init_angles = [0 for i in range(len(points))]
    thetas = numpy.array(init_angles, dtype=float)
    inlet_str= "-A*(x[1]-x_l)*(x[1]-x_u)"
    inlet_data = [[Expression((inlet_str, "0"), x_l=0.15, x_u=0.25,
                              A=250, degree=5), 1],
                  [Expression((inlet_str, "0"), x_l=0.73, x_u=0.83,
                              A=250, degree=5), 2]]

    pre = "meshes/"
    meshes = [pre+"multimesh_0.xdmf"] +  [pre+"multimesh_1.xdmf"]*len(points)
    mfs = [pre+"mf_0.xdmf"] + [pre+"mf_1.xdmf"]*len(points)
    solver = StokesSolver(points, thetas, meshes, mfs, inlet_data)

    # Save initial solution
    J_init = solver.eval_J(init_angles)
    print("Initial Theta")
    print(init_angles)
    print("Initial J: %.2e" % J_init)
    solver.save_state()

    
    # Initialize optimizer
    optimizer = IpoptAngle(len(thetas),solver.eval_J, solver.eval_dJ)
    optimizer.nlp.int_option('max_iter',20)
    #optimizer.nlp.num_option("tol", 1e-6)
    # optimizer.nlp.num_option("acceptable_tol", 1e-2)
    # optimizer.nlp.int_option("acceptable_iter", 3)
    optimizer.nlp.int_option("print_level",5)
    # optimizer.nlp.num_option("mu_init", 1)
    # optimizer.nlp.num_option("mu_min", 1)
    # optimizer.nlp.num_option("mu_max", 1)
    optimizer.nlp.str_option("output_file", "StokesOpt.txt") 
    # Solve optimization problem
    opt_theta = optimizer.solve(thetas)

    # # Save initial solution
    # J_init = solver.eval_J(init_angles)
    # print("Initial Theta")
    # print(init_angles)
    # print("Initial J: %.2e" % J_init)
    # solver.save_state()

    # Save final solution
    print("Optimal Theta")
    dJ_opt = solver.eval_dJ(opt_theta)
    J_opt = solver.eval_J(opt_theta)
    print("Optimal J: %.5f" % J_opt)
    print("Final gradient: ", dJ_opt)
    solver.save_state()
