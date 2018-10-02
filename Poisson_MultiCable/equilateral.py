import numpy
from IPython import embed
lmb_metal = 205.      # Heat coefficient aluminium
lmb_insulation = 0.03 # Heat coefficient of plastic
lmb_air = 0.15        # Heat coefficient of brick
c1 = numpy.array([0, 0.45])
c2 = numpy.array([-0.4, -0.15])
c3 = numpy.array([0.2,-0.4])
cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
scales = numpy.array([1,1,1])   
sources = numpy.array([5,15,5])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_air, sources)
opt = MultiCableOptimization(3, scales, MC.eval_J, MC.eval_dJ)

opt.nlp.num_option('obj_scaling_factor',1e-1) 
opt.nlp.int_option('max_iter', 100)
opt.nlp.num_option('acceptable_tol', 1e-2)
opt.nlp.num_option('tol', 1e-2)
opt_sol = opt.solve(cable_positions)

MC.eval_J(opt_sol)
from dolfin import plot
import matplotlib.pyplot as plt 
for i in range(MC.T.num_parts()):
    File("output/new_ipopt%d.pvd" % i) << MC.T.part(i)
