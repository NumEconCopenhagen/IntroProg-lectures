import numpy as np
from scipy import interpolate
from scipy import optimize

class ConsumptionSavingModel:

    def __init__(self, mp):
        ''' Initialize the model object

        Args:
            mp (SimpleNamespace) : model parameters
        
        Returns
            (ConsumptionSavingModel): model object
        '''
        
        # a. Parse parameters
        self.rho = mp.rho
        self.kappa = mp.kappa
        self.nu = mp.nu
        self.r0 = mp.r0
        self.beta = mp.beta
        self.Delta_y = mp.Delta_y
        self.y_prb_low = mp.y_prb_low
        self.Delta_r = mp.Delta_r
        self.n_r = mp.n_r 
        
        # b. Containers
        self.sim_m1 = []

    def utility(self, c):
        ''' Calculate flow utility of consumption level c

        Args:
            c (ndarray): level of consumtion

        Returns:
            (ndarray): flow utility of consumption
        '''

        return (c**(1-self.rho))/(1-self.rho)

    def bequest(self, m, c):
        ''' Calculate flow utility of leaving bequest given residual consumption

        Args:
            m (ndarray): cash-on-hand
            c (ndarray): level of consumtion

        Returns:
            (ndarray): utility of bequests
        '''

        return (self.nu*(m-c + self.kappa)**(1-self.rho))/(1-self.rho)

    def r_outcomes(self):
        ''' Create set of possible outcomes of interest rate and corresponding probabilites
        
        Returns:
            rs (ndarray): set of possible interest rate realizations
            r_prb (ndarray): probabilities corresponding to each interest rate realization
        '''

        # a. Create set of possible outcomes of r
        d1 = np.arange(1, self.n_r+1)
        d2 = np.flip(d1)*(-1)
        d = np.concatenate((d2, d1))
        rs = self.r0 + d*self.Delta_r

        # b. Uniform probability of each outcome of r
        r_prb = np.repeat(1.0 / (2*self.n_r), 2*self.n_r)

        return rs, r_prb

    def v2(self, c2, m2):
        ''' Compute state specific value of consumption choice and bequests in period 2

        Args:
            c2 (ndarray): level of consumtion in period 2
            m2 (ndarray): cash-on-hand in period 2

        Returns:
            (ndarray): value of comsumption and bequests
        '''

        return self.utility(c2) + self.bequest(m2,c2)

    def v1(self, c1, m1, v2_interp):
        ''' Compute state specific value of consumption choice in period 1

        Args:
            c1 (ndarray): level of consumtion in period 1
            m1 (ndarray): cash-on-hand in period 1
            v2_interp (RegularGridInterpolator): interpolator between m in period 2 and value function

        Returns:
            (ndarray): state specific value of consumption choice in period 1
        '''

        # a.1 Initialize variables
        expected_v2 = 0.0
        low_y = 1 - self.Delta_y
        high_y = 1 + self.Delta_y

        # a.2 Create set of possible outcomes of r 
        rs, r_prb = self.r_outcomes()

        # a.3 Assets at the end of period 1
        a1 = m1 - c1

        # b. Compute expectation of v2 given the set of possible interest rate and income realizations 
        for r,prb in zip(rs, r_prb):
            m2_low_y = (1+r)*a1 + low_y
            v2_low_y = self.y_prb_low*v2_interp([m2_low_y])

            m2_high_y = (1+r)*a1 + high_y
            v2_high_y = (1-self.y_prb_low)*v2_interp([m2_high_y])

            expected_v2 += prb*(v2_low_y + v2_high_y)

        # c. Return value v1 of consumption c1 and expected v2
        return self.utility(c1) + self.beta*expected_v2

    def solve_period_2(self):
        ''' Solve the consumption problem of period 2

        Returns:
            m2s (ndarray): cash-on-hand levels in period 2
            v2s (ndarray): value function in period 2
            c2s (ndarray): consumption function in period 2 (ie policy function)
        '''

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
            result = optimize.minimize(obj,[x0],method='L-BFGS-B',
            bounds=((1e-8,m2),))

            # iv. save
            v2s[i] = -result.fun
            c2s[i] = result.x
            
        return m2s,v2s,c2s
    
    def solve_period_1(self, v2_interp):
        ''' Solve the consumption problem of period 1

        Args:
            v2_interp (RegularGridInterpolator): interpolator between m in period 2 and value function

        Returns:
            m1s (ndarray): cash-on-hand levels in period 1
            v1s (ndarray): value function in period 1
            c1s (ndarray): consumption function in period 1 (ie policy function)
        '''

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
        ''' Solve the consumption savings problem over all periods

        Returns:
            m1 (ndarray): cash-on-hand levels in period 1
            v1 (ndarray): value function in period 1
            c1 (ndarray): optimal consumption function in period 1 (ie policy function)
            m2 (ndarray): cash-on-hand levels in period 2
            v2 (ndarray): value function in period 2
            c2 (ndarray): optimal consumption function in period 2 (ie policy function)
        '''

        # a. solve period 2
        m2, v2, c2 = self.solve_period_2()

        # b. construct interpolator
        v2_interp = interpolate.RegularGridInterpolator([m2], v2,
                                                        bounds_error=False, fill_value=None)

        # b. solve period 1
        m1, v1, c1 = self.solve_period_1(v2_interp)

        return m1, c1, v1, m2, c2, v2

    def simulate(self):
        ''' Simulate choices in period 1 and 2 based on model solution and random draws of income and interest rate. 
        
        Note: the parameters governing the random draws of income and interest rate are provided in the object mp when initializing the model.

        Returns:
            sim_c1 (ndarray): simulated consumption choices in period 1
            sim_c2 (ndarray): simulated consumption choices in period 2
        '''

        # a. solve the model at current parameters
        m1, c1, _, m2, c2, _ = self.solve()

        # b. construct interpolaters
        c1_interp = interpolate.RegularGridInterpolator([m1], c1,
                                                        bounds_error=False, fill_value=None)

        c2_interp = interpolate.RegularGridInterpolator([m2], c2,
                                                        bounds_error=False, fill_value=None)
    
        # c. sim period 1 based on draws of initial m and solution
        sim_c1 = c1_interp(self.sim_m1)
        sim_a1 = self.sim_m1-sim_c1

        # d. transition to period 2 m based on random draws of income and interest rate
        y2_down = 1-self.Delta_y
        y2_up = 1+self.Delta_y
        y2 = np.random.choice([y2_down, y2_up], 
                                p=[self.y_prb_low, (1-self.y_prb_low)], size=(sim_a1.shape))

        # e. Get the set of possible realizations of r for the current model and the corresponding probabilities.
        rs, r_prb = self.r_outcomes()

        # f. Use distribution of interest rate outcomes to create sample for simulation
        r = np.random.choice(list(rs), 
                                p=list(r_prb), size=(sim_a1.shape))

        # g. Based on random draws of income and interest rate, simulate period 2 cash on hand
        sim_m2 = (1+r)*sim_a1 + y2
            
        # h. Simulate period 2 consumption choice based on model solution 
        sim_c2 = c2_interp(sim_m2)

        return sim_c1, sim_c2, sim_m2

