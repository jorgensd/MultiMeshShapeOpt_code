from IPython import embed
import numpy
from scipy.optimize import minimize, NonlinearConstraint


class MultiCableOptimization():
    
    def __init__(self, num_cables, cable_scales):      
        self.g_scale = 1 # Scaling of coefficient for gradient constraint
        self.outer_radius = 1.2 # Radius of background cable
        self.inner_radius = 0.3 # Radius for each inner cable
        self.distance_from_outer = 0.05 # Distance for each cable to outer boundary
        self.distance_from_internal = 0.025 # Distance between each internal cable
        self.num_cables = num_cables
        self.inner_radius *= cable_scales
        self.max_radius = self.outer_radius - self.inner_radius\
                          - self.distance_from_outer
        self.sympy_g(self.num_cables)
        self.sympy_jac_g()
        self.constraint_init()

    def eval_g(self, x):
        """ Evaluate inequality constraint, g(x) <= 0, """
        return self.replace_sympy_g(x)


    def eval_jac_g(self, cable_positions, flag):
        """ The constraint Jacobian:
        flag = True  means 'tell me the sparsity pattern';
        flag = False means 'give me the Jacobian'.
        """
        if flag:
            # FIXME: Pass in a dense matrix
            # (could be optimised, but we don't care for this small problem).
            nvar = len(cable_positions)
            n_cables = int(nvar/2)
            ncon = int(n_cables + n_cables*(n_cables-1)/2)
            rows = []
            for i in range(ncon):
                rows += [i] * nvar
            cols = list(range(nvar)) * ncon
            return (numpy.array(rows), numpy.array(cols))
        else:
            return self.replace_sympy_jac_g(cable_positions)

    
    def sympy_g(self, num_cables):
        """ Creates the maximum distance and no-collision constraint 
        for each sub-cables with sympy, returning the g in inequality g<=0 and
        an array z of the control variables
        """
        x = list(sympy.symbols("x0:%d" % num_cables))
        y = list(sympy.symbols("y0:%d" % num_cables))
        g,z= [],[]
        for i in range(num_cables):
            g.append(x[i]**2 + y[i]**2 - self.max_radius[i]**2)
        for i in range(num_cables):
            z.append(x[i])
            z.append(y[i])
            for j in range(i+1,num_cables):
                xi, yi = x[i], y[i]
                xj, yj = x[j], y[j]
                int_radius = self.inner_radius[i] + self.inner_radius[j]\
                             + self.distance_from_internal
                # FIXME: Max_int_radius should change with each cable radius
                g.append(int_radius**2 - (xi - xj)**2 - (yi - yj)**2)
        self.g = sympy.Matrix(g)
        self.z = sympy.Matrix(z)
    

    def sympy_jac_g(self):
        """ Creates jacobian for constraints """
        self.jac_g = self.g.jacobian(self.z)


    def replace_sympy_g(self, positions):
        """ Create numpy array of constraint at current positions """
        g_numpy = self.g.subs([(self.z[i], positions[i])
                               for i in range(2*self.num_cables)])
        return self.g_scale*numpy.array(g_numpy).astype(float)


    def replace_sympy_jac_g(self, positions):
        """ Create numpy array of jacobian of constraint at
        current positions """    
        numpy_jac_g = self.jac_g.subs([(self.z[i], positions[i])
                                       for i in range(2*self.num_cables)])
        return self.g_scale*numpy.array(numpy_jac_g).astype(float)

    def constraint_init(self):
        nvar = int(2*self.num_cables)
        ncon = int(self.num_cables + self.num_cables*(self.num_cables-1)/2)
        low_var = -numpy.inf*numpy.ones(nvar,dtype=float)
        up_var = numpy.inf*numpy.ones(nvar, dtype=float)
        up_var[0], low_var[0] = 0, 0
        inf_con = numpy.inf*numpy.ones(ncon, dtype=float)
        zero_con = numpy.zeros(ncon, dtype=float)
        self.non_lin_constraint = NonlinearConstraint(self.eval_g,
                                                      low_var,
                                                      up_var,
                                                      jac=self.eval_jac_g,
                                                      keep_feasible=True)

    def solve(self, cable_positions, eval_J, eval_dJ, callback,
              options={}):
        result = minimize(eval_J, thetas,
                          jac=eval_dJ,
                          method = 'Newton-CG',
                          callback=callback,
                          options={"maxiter":1,"xtol":1e-2, "disp":True})
        
        return result

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


if __name__ == "__main__":
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

    MC = MultiCable(scales, cable_positions, lmb_metal, lmb_insulation,
                      lmb_air, sources)
 
    
    # Save initial solution
    J_init = MC.eval_J(cable_positions)
    print("Initial J: %.2e" % J_init)
    MC.save_state()

    # Optimize
    optimizer = MultiCableOptimization(MC.num_cables, scales)
    result = optimizer.solve(cable_positions, MC.eval_J, MC.eval_dJ,
                             MC.callback)
    
    # Save final solution
    print("Optimal positions")
    print(result["x"])
    dJ_opt = MC.eval_dJ(result["x"])
    J_opt = MC.J
    print("Optimal J: %.5f" % J_opt)
    print("Final gradient: ", dJ_opt)
    MC.save_state()
