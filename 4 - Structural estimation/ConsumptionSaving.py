import numpy as np
from scipy import interpolate
from scipy import linalg
from scipy import optimize

class ConsumptionSavingModel:

    def __init__(self, par):
        self.par = par
        self.sim_m1 = []
        self.data_m1 = []
        pass

    def utility(self,c):
        return c**(1-self.par.rho)/(1-self.par.rho)

    def bequest(self,m,c):
        return self.par.nu*(m-c+self.par.kappa)**(1-self.par.rho)/(1-self.par.rho)

    def v2(self,c2,m2):
        return self.utility(c2) + self.bequest(m2,c2)

    def v1(self,c1,m1,v2_interp):
        
        # a. v2 value, if low income
        m2_low = (1+self.par.r)*(m1-c1) + 1-self.par.Delta
        v2_low = v2_interp([m2_low])[0]
        
        # b. v2 value, if high income
        m2_high = (1+self.par.r)*(m1-c1) + 1+self.par.Delta
        v2_high = v2_interp([m2_high])[0]
        
        # c. expected v2 value
        prob_low = 0.5
        prob_high = 0.5
        expected_v2 = prob_low*v2_low + prob_high*v2_high
        
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
        m1s = np.linspace(1e-8, 4, 100)
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

        # a. solve period 2
        m2, v2, c2 = self.solve_period_2()

        # b. construct interpolator
        v2_interp = interpolate.RegularGridInterpolator([m2], v2,
                                                        bounds_error=False, fill_value=None)

        # b. solve period 1
        m1, v1, c1 = self.solve_period_1(v2_interp)

        return m1, c1, m2, c2

    def simulate(self):

        # a. solve the model at current parameters
        m1, c1, m2, c2 = self.solve()

        # b. construct interpolaters
        c1_interp = interpolate.RegularGridInterpolator([m1], c1,
                                                        bounds_error=False, fill_value=None)

        c2_interp = interpolate.RegularGridInterpolator([m2], c2,
                                                        bounds_error=False, fill_value=None)
    
        # c. sim period 1 based on draws of initial m and solution
        sim_c1 = c1_interp(self.sim_m1)
        sim_a1 = self.sim_m1-sim_c1

        # d. transition to period 2 m based on random draws
        sim_m2 = (1+self.par.r)*sim_a1 + np.random.choice([0.5, 1.5], p=[0.5, 0.5], size=(sim_a1.shape))

        # e. sim period 2 consumption choice based on model solution and sim_m2
        sim_c2 = c2_interp(sim_m2)

        return sim_c1, sim_c2

