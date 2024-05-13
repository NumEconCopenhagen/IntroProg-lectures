# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})
from scipy import interpolate


# Define the model
class ConsumerModelClass():

    def __init__(self, **kwargs):
        '''
        Initialize the model with default parameters
        kwargs allow any parameter in the par namespace to be overridden by the user
        '''

        self.par = par = SimpleNamespace() # Create a namespace object for parameters
        self.sol = sol = SimpleNamespace() # Create a namespace object for solution results
        self.sim = sim = SimpleNamespace() # Create a namespace object for simulation results (not always neccesary)

        # Set default parameters
        self.setup()

        # Update parameters with user input
        for key, value in kwargs.items():
            setattr(par, key, value)

        # Allocate arrays simulation (if needed)
        self.allocate()


    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par

        # Model parameters
        par.rho = 8
        par.kappa = 0.5
        par.nu = 0.1
        par.r = 0.04
        par.beta = 0.94
        par.Delta = 0.5

        
        # Helpful variables 
        par.approx0 = 1e-8

        par.vecN = 500



    def allocate(self):
        '''
        Allocate arrays for simulation
        '''

        par = self.par
        sol = self.sol

        # solution vecs
        vecs = ['m1','m2','c1','c2','v1','v2']
        for vec in vecs:
            setattr(sol,f'{vec}_vec',np.empty(par.vecN))
    
        sol.m2_vec[:] = np.linspace(par.approx0,5,par.vecN)
        sol.m1_vec[:] = np.linspace(par.approx0,4,par.vecN)
            
    def utility(self, c):
        par = self.par
        return c**(1-par.rho)/(1-par.rho)

    def bequest(self, m, c):
        par = self.par
        return par.nu*(m-c+par.kappa)**(1-par.rho)/(1-par.rho)

    def v2(self, c2, m2):
        return self.utility(c2) + self.bequest(m2, c2)

    def v1(self, c1, m1):
        par = self.par
        sol = self.sol 

        # a. v2 value, if low income
        m2_low = (1+par.r)*(m1-c1) + 1-par.Delta
        v2_low = sol.v2_interp([m2_low])[0]
        
        # b. v2 value, if high income
        m2_high = (1+par.r)*(m1-c1) + 1+par.Delta
        v2_high = sol.v2_interp([m2_high])[0]
        
        # c. expected v2 value
        v2 = 0.5*v2_low + 0.5*v2_high
        
        # d. total value
        return self.utility(c1) + par.beta*v2
    
    def v1_alt(self,c1,m1):
        par = self.par 
        sol = self.sol 

        probs = np.array([0.1,0.4,0.4,0.1])
        y2s = np.array([1-np.sqrt(par.Delta),1-par.Delta,1+par.Delta,1+np.sqrt(par.Delta)])

        # 
        m2_vec = (1+par.r)*(m1-c1) +y2s 
        v2s =  sol.v2_interp(m2_vec)
        
        # expected v2 value
        v2 = probs@v2s
        
        # total value
        return self.utility(c1) + par.beta*v2


    def solve_period_2(self):
        par = self.par
        sol = self.sol        

        # b. solve for each m2 in grid
        for i,m2 in enumerate(sol.m2_vec):

            # i. objective
            obj = lambda x: -self.v2(x[0],m2)

            # ii. initial value (consume half)
            x0 = m2/2

            # iii. optimizer
            result = optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((par.approx0/100,m2),))

            # iv. save
            sol.v2_vec[i] = -result.fun
            sol.c2_vec[i] = result.x[0]
            
        

    def solve_period_1(self,v1_type='original'):
        par = self.par
        sol = self.sol

        # b. solve for each m1 in grid
        for i,m1 in enumerate(sol.m1_vec):
            
            # i. objective
            if v1_type == 'original':
                obj = lambda x: -self.v1(x[0],m1)
            elif v1_type == 'alt':
                obj = lambda x: -self.v1_alt(x[0],m1)
            
            # ii. initial guess (consume half)
            x0 = m1/2
            
            # iii. optimize
            result = optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((par.approx0/100,m1),))
            
            # iv. save
            sol.v1_vec[i] = -result.fun
            sol.c1_vec[i] = result.x[0]
    

    def solve(self,v1_type='original'):
        par = self.par
        sol = self.sol

        # a. solve period 2
        self.solve_period_2()

        # b. construct interpolator
        sol.v2_interp = interpolate.RegularGridInterpolator((sol.m2_vec,), sol.v2_vec,
                                                    bounds_error=False,fill_value=None)
        # c. solve period 1
        self.solve_period_1(v1_type=v1_type)


    def plot_c1(self,add_alt=False):
        par = self.par
        sol = self.sol

        fig,ax  = plt.subplots()
        ax.plot(sol.m1_vec,sol.c1_vec,label='Original')
        ax.set_xlabel('$m_1$')
        ax.set_ylabel('$c_1$')
        ax.set_title('Consumption in period 1')

        if add_alt:
            # v2_interp is the same for this new utility function in period 1
            self.solve_period_1(v1_type='alt')
            ax.plot(sol.m1_vec,sol.c1_vec,label='Alternative utility function')
            ax.legend()
        plt.show()