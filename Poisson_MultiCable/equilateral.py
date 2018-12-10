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
sources = numpy.array([10,10,10])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_air, sources)
from dolfin import plot, File
num_cells = 0 
for i in range(MC.multimesh.num_parts()):
    num_cells +=(len(MC.multimesh.cut_cells(i))
                 +len(MC.multimesh.uncut_cells(i)))
print("#Cells: %d" %num_cells)
MC.eval_J(cable_positions)
print("initial J %.3e" % MC.J)


opt = MultiCableOptimization(3, scales, MC.eval_J, MC.eval_dJ)
opt.nlp.int_option('max_iter', 30)

opt_sol = opt.solve(cable_positions)


compute_angles(opt_sol)
MC.eval_J(opt_sol)
print("Optimal J: %.2e" % MC.J)
def plot_init_opt(init, opt,filename,colorbar=True):
    plt.subplot(1,2,1)
    MC.eval_J(init)
    MC.T.vector().get_local().copy()
    
    T_max_arr = MC.T.vector().get_local()
    inactive = MC.T.function_space().dofmap().inactive_dofs(MC.multimesh, 0)
    T_stripped = [val for idx, val in enumerate(T_max_arr)
                  if idx not in inactive]
    T_max = numpy.max(T_stripped)
    MC.eval_J(opt)
    T_min_arr = MC.T.vector().get_local()   
    inactive = MC.T.function_space().dofmap().inactive_dofs(MC.multimesh, 0)
    T_stripped = [val for idx, val in enumerate(T_min_arr)
                  if idx not in inactive]
    T_min = numpy.min(T_stripped)
    MC.eval_J(init)
    
    # Strip T of inactive dofs values
    N = 512 # Number of levels used in contour plot
    
    plot_cable(init,N,T_min,T_max)

    plt.subplot(1,2,2)
    MC.eval_J(opt)
    
    plot_cable(opt,N, T_min, T_max)
    m = plt.cm.ScalarMappable()
    m.set_array(T_stripped)
    fig = plt.gcf()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.81, 0.32, 0.02, 0.34])
    cbar = plt.colorbar(m,orientation="vertical",
                        ticks=numpy.linspace(T_min,T_max, 4),
                        boundaries=numpy.linspace(T_min,T_max,N),pad=0.01,
                        format="%.2f", cax=cbar_ax)

    cbar.set_clim(T_min,T_max)
    plt.savefig(filename, pad_inches=0, dpi=250,
                bboxinches="tight")
    import os
    os.system("convert " + filename + " -trim " + filename)
    
def plot_cable(positions,N, T_min, T_max):
    fig = plt.gcf()
    ax = fig.gca()
    norm = mpl.colors.Normalize(vmin=T_min, vmax=T_max)
    
    p0 = plot(MC.T.part(0),norm=norm,zorder=0,
              levels=numpy.linspace(T_min,T_max,N))

    for i in range(1,MC.num_cables+1):
        pi = plot(MC.T.part(i), zorder=i,norm=norm)
        rubber_radius = 0.85*opt.inner_radius[i-1]
        metal_radius = 2./3*opt.inner_radius[i-1]
        iso_circle = plt.Circle((positions[2*(i-1)], positions[2*(i-1)+1]),
                            radius=rubber_radius, color="k",
                                fill=False,zorder=10+i,linewidth=0.5)
        metal_circle = plt.Circle((positions[2*(i-1)], positions[2*(i-1)+1]),
                                  radius=metal_radius, color="k",
                                  fill=False,zorder=10+i,linewidth=0.5)

        ax.add_patch(iso_circle)
        ax.add_patch(metal_circle)
    plt.axis("off")
    return T_min, T_max
    
plot_init_opt(cable_positions, opt_sol,"output/equilateral.png")

for i in range(MC.multimesh.num_parts()):
    File("single_mesh/output/opt_mc%d.pvd" %i) << MC.multimesh.part(i)
list_timings(TimingClear.clear, [TimingType.wall])
