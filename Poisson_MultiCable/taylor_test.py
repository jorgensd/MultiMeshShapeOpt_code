import numpy
import matplotlib.pyplot as plt
from IPython import embed

def convergence_rates(E_values, eps_values):
    r = []
    for i in range(1, len(eps_values)):
        r.append(numpy.log(E_values[i]/E_values[i-1])/
                 numpy.log(eps_values[i]/eps_values[i-1]))
    return r

lmb_metal = 205.   # Heat coefficient aluminium
lmb_insulation = 0.33 # Heat coefficient of plastic
lmb_air = 0.33   # Heat coefficient of brick
c1 = numpy.array([0, 0.2,-0.2,-0.4])
scales = numpy.array([1,1])   
sources = numpy.array([25,25])
from MultiCable import *

MC = MultiCable(scales, c1, lmb_metal, lmb_insulation,
                      lmb_air, sources)
num_cells = 0 
for i in range(MC.multimesh.num_parts()):
    num_cells += MC.multimesh.part(i).num_cells()

from dolfin import plot, File

perturbation= numpy.array([0,0.65,-0.1,0.2])
dJ = MC.eval_dJ(c1)
dJp = numpy.dot(dJ, perturbation)
J = MC.J

epsilon = [0.5**(i) for i in range(10)]
res_0 = []
res_1 = []
outputs = [File("output/taylor%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]
for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)
for eps in epsilon:
    J_eps = MC.eval_J(c1+eps*perturbation)
    res_0.append(numpy.abs(J_eps-J))
    res_1.append(numpy.abs(J_eps-J-eps*dJp))
    for i in range(MC.multimesh.num_parts()):
        outputs[i] << MC.T.part(i)

print("#Cells", num_cells)
print(' '.join('{:1.5e}'.format(k) for k in epsilon))
print(' '.join('{:1.5e}'.format(k) for k in res_0))
print(' '.join('{:1.5e}'.format(k) for k in res_1))
print(' '.join('{:1.5e}'.format(k) for k in convergence_rates(res_0,epsilon)))
print(' '.join('{:1.5e}'.format(k) for k in convergence_rates(res_1,epsilon)))

