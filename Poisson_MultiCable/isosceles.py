import numpy

def compute_angles(cable_positions):
    """
    Compute angles between three cables
    """
    c1 = cable_positions[0:2]
    c2 = cable_positions[2:4]
    c3 = cable_positions[4:6]
    a1 = (numpy.arccos(numpy.dot(c1-c2,c1-c3)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c1-c3,c1-c3))))/(2*numpy.pi)*360)
    a2 = (numpy.arccos(numpy.dot(c1-c2,c3-c2)/
                       (numpy.sqrt(numpy.dot(c1-c2,c1-c2)
                                   *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)
    a3 = (numpy.arccos(numpy.dot(c3-c1,c3-c2)/
                      (numpy.sqrt(numpy.dot(c3-c1,c3-c1)
                                  *numpy.dot(c3-c2,c3-c2))))/(2*numpy.pi)*360)
    print("Angles between the three cables: %.2f, %.2f, %.2f" %(a1,a2,a3))

lmb_metal = 205.      # Heat coefficient aluminium
lmb_insulation = 0.03 # Heat coefficient of plastic
lmb_air = 0.15        # Heat coefficient of brick
c1 = numpy.array([-0.2, 0.45])
c2 = numpy.array([-0.4, -0.3])
c3 = numpy.array([0.3,-0.1])
cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
compute_angles(cable_positions)

scales = numpy.array([1,1,1])   
sources = numpy.array([10,5,5])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_air, sources)
from dolfin import plot, File
outputs = [File("output/isosceles%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]
MC.eval_J(cable_positions)
for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)


opt = MultiCableOptimization(3, scales, MC.eval_J, MC.eval_dJ)
# opt.nlp.num_option('obj_scaling_factor',1e-2) 
# opt.nlp.int_option('max_iter', 50)
# opt.nlp.num_option('acceptable_tol', 1e-3)
# opt.nlp.num_option('tol', 5e-4)

opt_sol = opt.solve(cable_positions)
compute_angles(opt_sol)
MC.eval_J(opt_sol)

for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)
