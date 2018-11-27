import numpy
import matplotlib.pyplot as plt
from IPython import embed

def convergence_rates(E_values, eps_values):
    r = []
    for i in range(1, len(eps_values)):
        r.append(numpy.log(E_values[i]/E_values[i-1])/
                 numpy.log(eps_values[i]/eps_values[i-1]))
    return r


lmb_metal = 210.      # Heat coefficient aluminium
lmb_insulation = 0.44 # Heat coefficient of plastic
lmb_air = 0.03   # Heat coefficient of brick
c1 = numpy.array([0.02, 0.1])


scales = numpy.array([1])   
sources = numpy.array([10])

from MultiCable import *

MC = MultiCable(scales, c1, lmb_metal, lmb_insulation,
                      lmb_air, sources)

from dolfin import plot, File
outputs = [File("output/taylor%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]

perturbation= numpy.array([0.5,0.5])
dJ = MC.eval_dJ(c1)
dJp = numpy.dot(dJ, perturbation)

J = MC.J
for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)

epsilon = [0.5**(i) for i in range(5)]
res_0 = []
res_1 = []
for eps in epsilon:
    J_eps = MC.eval_J(c1+eps*perturbation)
    for i in range(MC.multimesh.num_parts()):
        outputs[i] << MC.T.part(i)
    res_0.append(numpy.abs(J_eps-J))
    res_1.append(numpy.abs(J_eps-J-eps*dJp))
print("FD residuals")
print(["%.2e" % i for i in res_0])
print("FD rate")
print(["%.2f" % i for i in convergence_rates(res_0, epsilon)])
print("1st order residuals")
print(["%.2e" % i for i in res_1])
print("1st order rate")
print(["%.2f" % i for i in convergence_rates(res_1, epsilon)])
