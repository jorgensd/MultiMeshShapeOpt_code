import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython import embed

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
lmb_air = 0.33        # Heat coefficient of brick
c1 = numpy.array([0, 0.45])
c2 = numpy.array([-0.4, -0.15])
c3 = numpy.array([0.2,-0.4])
cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1]])
compute_angles(cable_positions)

scales = numpy.array([1,1,1])   
sources = numpy.array([1,1,1])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_air, sources)
from dolfin import plot, File
outputs = [File("output/equilateral%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]
MC.eval_J(cable_positions)

for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)


opt = MultiCableOptimization(3, scales, MC.eval_J, MC.eval_dJ)
opt.nlp.int_option('max_iter', 30)
opt.nlp.num_option('tol', 1e-6)


def plot_cable(positions):
    MC.eval_J(positions)
    fig = plt.figure()
    ax = fig.gca()
    # Strip T of inactive dofs values
    T_arr = MC.T.vector().get_local()
    inactive = MC.T.function_space().dofmap().inactive_dofs(MC.multimesh, 0)
    T_stripped = [val for idx, val in enumerate(T_arr)
                       if idx not in inactive]
    T_min = numpy.min(T_stripped)
    T_max = numpy.max(T_stripped)
    norm = mpl.colors.Normalize(vmin=T_min, vmax=T_max)
    N = 512 # Number of levels used in contour plot
    p0 = plot(MC.T.part(0),norm=norm,zorder=0,
              levels=numpy.linspace(T_min,T_max,N))

    for i in range(1,MC.num_cables+1):
        pi = plot(MC.T.part(i), zorder=i,norm=norm)
        rubber_radius = 0.85*opt.inner_radius[i-1]
        metal_radius = 2./3*opt.inner_radius[i-1]
        iso_circle = plt.Circle((positions[2*(i-1)], positions[2*(i-1)+1]),
                            radius=rubber_radius, color="k",
                            fill=False,zorder=10+i)
        metal_circle = plt.Circle((positions[2*(i-1)], positions[2*(i-1)+1]),
                            radius=metal_radius, color="k",
                            fill=False,zorder=10+i)

        ax.add_patch(iso_circle)
        ax.add_patch(metal_circle)
    ax = fig.gca()
    # Creating scaled colorbar from stripped data
    m = plt.cm.ScalarMappable()
    m.set_array(T_stripped)
    cbar = plt.colorbar(m,boundaries=numpy.linspace(T_min,T_max,N))
    cbar.set_clim(T_min,T_max)
    plt.savefig("output/init_cable.png")
    exit(1)

plot_cable(cable_positions)
exit(1)

opt_sol = opt.solve(cable_positions)
compute_angles(opt_sol)
MC.eval_J(opt_sol)

# for i in range(MC.multimesh.num_parts()):
#     outputs[i] << MC.T.part(i)

