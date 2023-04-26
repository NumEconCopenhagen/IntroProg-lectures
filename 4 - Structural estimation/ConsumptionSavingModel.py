from types import SimpleNamespace

import numpy as np
from scipy import interpolate
from scipy import optimize

class ConsumptionSavingModelClass:

    def __init__(self):

        # a. namespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # b. parameters
        par.rho = 2.0
        par.kappa = 0.5
        par.nu = 10.0
        par.r = 0.04
        par.beta = 0.94        
        par.simN = 1_000

    def utility(self,c):
        return c**(1-self.par.rho)/(1-self.par.rho)

    def bequest(self,m,c):
        return self.par.nu*(m-c+self.par.kappa)**(1-self.par.rho)/(1-self.par.rho)

    def v2(self,c2,m2):
        return self.utility(c2) + self.bequest(m2,c2)

    def v1(self,c1,m1,v2_interp):
        
        # a. v2 value, if low income
        m2_low = (1+self.par.r)*(m1-c1) + 0.5
        v2_low = v2_interp([m2_low])[0]
        
        # b. v2 value, if high income
        m2_high = (1+self.par.r)*(m1-c1) + 1.5
        v2_high = v2_interp([m2_high])[0]
        
        # c. expected v2 value
        expected_v2 = 0.5*v2_low + 0.5*v2_high
        
        # d. total value
        return self.utility(c1) + self.par.beta*expected_v2

    def solve_period_2(self):

        # a. grids
        m2s = np.linspace(1e-4,5,500)
        v2s = np.empty(500)
        c2s = np.empty(500)

        # b. solve for each m2 in grid
        for i,m2 in enumerate(m2s):

            # i. objective
            obj = lambda x: -self.v2(x[0],m2)

            # ii. initial value (consume half)
            x0 = m2/2

            # iii. optimizer
            result = optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((1e-8,m2),))

            # iv. save
            v2s[i] = -result.fun
            c2s[i] = result.x
            
        return m2s,v2s,c2s
    
    def solve_period_1(self, v2_interp):

         # a. grids
        m1s = np.linspace(1e-8,4,100)
        v1s = np.empty(100)
        c1s = np.empty(100)

        # b. solve for each m1s in grid
        for i, m1 in enumerate(m1s):

            # i. objective
            def obj(x): return -self.v1(x[0], m1, v2_interp)

            # ii. initial guess (consume half)
            x0 = m1/2

            # iii. optimize
            result = optimize.minimize(
                obj, [x0], method='L-BFGS-B', bounds=((1e-12, m1),))

            # iv. save
            v1s[i] = -result.fun
            c1s[i] = result.x[0]

        return m1s, v1s, c1s
    
    def solve(self):
        """ solve the model """

        sol = self.sol

        # a. solve period 2
        sol.m2, sol.v2, sol.c2 = self.solve_period_2()

        # b. construct interpolator
        v2_interp = interpolate.RegularGridInterpolator([sol.m2], sol.v2,
                                                        bounds_error=False, fill_value=None)

        # b. solve period 1
        sol.m1,sol.v1,sol.c1 = self.solve_period_1(v2_interp)

    def draw_random_values(self):
        """ draw random values for simulation """

        par = self.par
        sim = self.sim

        sim.m1 = np.fmax(np.random.normal(1,0.1,size=par.simN),0)
        sim.y2 = np.random.choice([0.5,1.5],p=[0.5,0.5],size=(par.simN))

    def simulate(self):
        """ simulate the model """

        par = self.par
        sol = self.sol
        sim = self.sim

        # a. construct interpolaters
        c1_interp = interpolate.RegularGridInterpolator([sol.m1], sol.c1,
                                                        bounds_error=False, fill_value=None)

        c2_interp = interpolate.RegularGridInterpolator([sol.m2], sol.c2,
                                                        bounds_error=False, fill_value=None)

        # b. sim period 1 based on draws of initial m and solution
        sim.c1 = c1_interp(sim.m1)
        sim.a1 = sim.m1-sim.c1

        # c. transition to period 2 m based on random draws
        sim.m2 = (1+par.r)*sim.a1 + sim.y2

        # d. sim period 2 consumption choice based on model solution and sim.m2
        sim.c2 = c2_interp(sim.m2)