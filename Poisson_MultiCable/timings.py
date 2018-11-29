import numpy
import matplotlib.pyplot as plt
from IPython import embed

lmb_metal = 205.   # Heat coefficient aluminium
lmb_insulation = 0.2 # Heat coefficient of plastic
lmb_air = 0.33   # Heat coefficient of brick
c1 = numpy.array([0, 0.1])


scales = numpy.array([1])   
sources = numpy.array([25])

from MultiCable import *

MC = MultiCable(scales, c1, lmb_metal, lmb_insulation,
                      lmb_air, sources)

outputs = [File("output/taylor%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]

dJ = MC.eval_dJ(c1)

perturbation= numpy.array([0,0.75])
epsilon = [0.5**(i) for i in range(12)]
for eps in epsilon:
    dJ_eps = MC.eval_dJ(c1+eps*perturbation)

list_timings(TimingClear.clear, [TimingType.wall])
embed()
dJp = numpy.dot(dJ, perturbation)
print(MC.J,dJ)
J = MC.J
