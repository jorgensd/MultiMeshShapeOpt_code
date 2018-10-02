import numpy

lmb_metal = [205.,205,205,205,205]
lmb_insulation = [0.06,0.03,0.02,0.02,0.02]
lmb_fill = 0.33 
c1 = numpy.array([-0.2, -0.4])
c2 = numpy.array([-0.4, 0.2])
c3 = numpy.array([0, 0.71])
c4 = numpy.array([0.52, 0.3])
c5 = numpy.array([0.45, -0.45])
cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1],
                               c4[0],c4[1],c5[0],c5[1]])

scales = numpy.array([0.9,1,0.6,1,0.8])   
sources = numpy.array([5,4,1,3,2])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_fill, sources)
from dolfin import plot, File
outputs = [File("output/fivecables%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]
MC.eval_J(cable_positions)
for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)


opt = MultiCableOptimization(5, scales, MC.eval_J, MC.eval_dJ)
# opt.nlp.num_option('obj_scaling_factor',1e-1) 
# opt.nlp.int_option('max_iter', 20)
# opt.nlp.num_option('acceptable_tol', 1e-3)
# opt.nlp.num_option('tol', 2.5e-4)

opt_sol = opt.solve(cable_positions)

MC.eval_J(opt_sol)

for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)
