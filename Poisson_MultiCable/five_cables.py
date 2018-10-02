import numpy

lmb_metal = [205.,205,205,205,205]
lmb_insulation = [0.03,0.05,0.03,0.04,0.06]
lmb_fill = 0.15 
c1 = numpy.array([-0.2, -0.4])
c2 = numpy.array([-0.4, 0.2])
c3 = numpy.array([0, 0.6])
c4 = numpy.array([0.36, 0.3])
c5 = numpy.array([0.45, -0.45])
cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1],
                               c4[0],c4[1],c5[0],c5[1]])

scales = numpy.array([0.5,0.75,0.6,0.9,1])   
sources = numpy.array([15,7,5,3,2])

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
opt.nlp.num_option('obj_scaling_factor',1e-2) 
opt.nlp.int_option('max_iter', 50)
# opt.nlp.num_option('acceptable_tol', 1e-3)
opt.nlp.num_option('tol', 5e-4)

opt_sol = opt.solve(cable_positions)

MC.eval_J(opt_sol)

for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)
