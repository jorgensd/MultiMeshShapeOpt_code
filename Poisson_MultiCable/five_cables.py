import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
lmb_metal = [205,205,205,205,205]
lmb_insulation = [0.03,0.12,0.06,0.04,0.02]
lmb_fill = 0.33
c1 = numpy.array([0, 0.6])
c2 = numpy.array([-0.4, 0.2])
c3 = numpy.array([-0.1, -0.4])
c4 = numpy.array([0.6, 0.4])
c5 = numpy.array([0.45, -0.45])

cable_positions = numpy.array([c1[0],c1[1],c2[0],c2[1],c3[0],c3[1],
                               c4[0],c4[1],c5[0],c5[1]])

scales = numpy.array([1,0.75,0.9,1,0.8])   
sources = numpy.array([10,5,2.5,5,10])

from MultiCable import *
from IpoptMultiCableSolver import *

MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_fill, sources)
from dolfin import plot, File
outputs = [File("output/fivecables%d.pvd" %i)
           for i in range(MC.multimesh.num_parts())]
MC.eval_J(cable_positions)
num_cells = 0 
for i in range(MC.multimesh.num_parts()):
    num_cells +=(len(MC.multimesh.cut_cells(i))
                 +len(MC.multimesh.uncut_cells(i)))
print("#Cells: %d" %num_cells)

for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)

n_cables = MC.multimesh.num_parts()-1
opt = MultiCableOptimization(n_cables, scales, MC.eval_J, MC.eval_dJ)
opt.nlp.int_option('max_iter', 50)
# opt.nlp.num_option('tol', 1e-8)
opt_sol = opt.solve(cable_positions)

MC.eval_J(opt_sol)
for i in range(n_cables):
    print("%.8f, %.8f" %(opt_sol[2*i], opt_sol[2*i+1]))
        
for i in range(MC.multimesh.num_parts()):
    outputs[i] << MC.T.part(i)

def plot_init_opt(init, opt,filename):
    plt.subplot(1,2,1)
    MC.eval_J(init)
    # Strip T of inactive dofs values
    N = 512 # Number of levels used in contour plot
    
    plot_cable(init,N)

    plt.subplot(1,2,2)
    MC.eval_J(opt)
    
    plot_cable(opt,N)
    plt.savefig(filename, pad_inches=0, dpi=250,
                bboxinches="tight")
    import os
    os.system("convert " + filename + " -trim " + filename)
    
def plot_cable(positions,N):
    fig = plt.gcf()
    ax = fig.gca()
    T_arr = MC.T.vector().get_local()
    inactive = MC.T.function_space().dofmap().inactive_dofs(MC.multimesh, 0)
    T_stripped = [val for idx, val in enumerate(T_arr)
                       if idx not in inactive]
    T_min = numpy.min(T_stripped)
    T_max = numpy.max(T_stripped)
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
    m = plt.cm.ScalarMappable()
    m.set_array(T_stripped)
    cbar = plt.colorbar(m,orientation="horizontal",
                        ticks=numpy.linspace(T_min,T_max, 4),
                        boundaries=numpy.linspace(T_min,T_max,N),pad=0.01,
                        format="%.2f")

    cbar.set_clim(T_min,T_max)
    plt.axis("off")

    
plot_init_opt(cable_positions, opt_sol,"output/five_cables.png")
