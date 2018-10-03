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
        low_var = -360*numpy.ones(nvar,dtype=float)
        up_var = 360*numpy.ones(nvar, dtype=float)


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
        self.nlp.num_option('bound_relax_factor', 0)

    def solve(self, angle):
        return self.nlp.solve(angle)[0]


from Poisson_solver import *
p = Point(1.25,0.875)
theta = numpy.array([-45], dtype=float)
m_names = ["meshes/multimesh_%d.xdmf" %i for i in range(2)]
f_names = ["meshes/mf_%d.xdmf" %i for i in range(2)]
fexp = Expression('10*x[0]*sin(x[0])*cos(x[1])', degree=4)
solver = PoissonSolver(p, theta[0], m_names, f_names, fexp)
File("output/firstmesh.pvd") << solver.T.part(1)
optimizer = IpoptAngle(1,solver.eval_J, solver.eval_dJ)
optimizer.nlp.int_option('max_iter',10)
optimizer.nlp.num_option('obj_scaling_factor',1e0)
opt_theta = optimizer.solve(theta)[0]

print(opt_theta)
File("output/tmp0Opt.pvd") << solver.T.part(0)
File("output/tmpOpt.pvd") << solver.T.part(1)

