from IPython import embed
from dolfin import *
from StokesSolver import *
import numpy

if __name__ == "__main__":
    Points = [0.25,0.25,0.5,0.25,0.75,0.25,
              0.25,0.5, 0.5,0.5, 0.75,0.5,
              0.25,0.75,0.5,0.75,0.75,0.75]
    points = [Point(Points[2*i],Points[2*i+1])
              for i in range(int(len(Points)/2))]
    # points = [Point(0.24,0.5)]
    init_angles = [0 for i in range(len(points))]
    thetas = numpy.array(init_angles, dtype=float)
    inlet_str= "-A*(x[1]-x_l)*(x[1]-x_u)"
    inlet_data = [[Expression((inlet_str, "0"), x_l=0.15, x_u=0.25,
                              A=250, degree=5), 1],
                  [Expression((inlet_str, "0"), x_l=0.73, x_u=0.83,
                              A=250, degree=5), 2]]

    pre = "meshes/"
    meshes = [pre+"multimesh_0.xdmf"] +  [pre+"multimesh_1.xdmf"]*len(points)
    mfs = [pre+"mf_0.xdmf"] + [pre+"mf_1.xdmf"]*len(points)
    solver = StokesSolver(points, thetas, meshes, mfs, inlet_data)

    # Save initial solution
    J_init = solver.eval_J(init_angles)
    print("Initial Theta")
    print(init_angles)
    print("Initial J: %.2e" % J_init)
    solver.save_state()

    # Initialize optimizer
    from scipy.optimize import minimize
    result = minimize(solver.eval_J, thetas,
                         jac=solver.eval_dJ,
                         method = 'Newton-CG',
                         callback=solver.callback,
                         options={"maxiter":50,"xtol":1e-2, "disp":True})

    # Save final solution
    print("Optimal Theta")
    print(result["x"])
    dJ_opt = solver.eval_dJ(result["x"])
    J_opt = solver.J
    print("Optimal J: %.5f" % J_opt)
    print("Final gradient: ", dJ_opt)
    solver.save_state()
