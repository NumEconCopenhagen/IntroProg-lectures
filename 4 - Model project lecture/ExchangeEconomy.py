
# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})


class ExchangeEconomyModel:
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

        # Simulate draws
        self.simulate()

        
    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par

        # a. parameters
        par.N = 50
        par.mu = np.array([3,2,1])
        par.Sigma = np.array([[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]])
        par.gamma = 0.8
        par.zeta = 1

        par.epsilon = 1e-4
        par.kappa = 4
       



    def allocate(self):
        '''
        Allocate arrays for simulation
        '''


        par = self.par
        sim = self.sim

        simvarnames = ['e1','e2','e3'] # Which arrays should be available in the simulation namespace

        for varname in simvarnames:
            sim.__dict__[varname] =  np.nan*np.ones(par.N) # Allocate the size of the arrays
        

        N3_vars = ['alphas','betas']
        for varname in N3_vars:
            sim.__dict__[varname] =  np.nan*np.ones((par.N,3)) 

    
    def simulate(self,seed=1234):
        '''
        Draw random shocks for your model

        It's a good idea to draw all the random shocks you want in one function and store them 
        You can then be explicit about when you change the seed and take a new draw.
        '''
        par = self.par
        sim = self.sim 
        
        np.random.seed(seed)

        # preferences
        sim.alphas[:] = np.exp(np.random.multivariate_normal(par.mu, par.Sigma, size=par.N))
        sim.betas[:] = sim.alphas/np.reshape(np.sum(sim.alphas,axis=1),(par.N,1))

        # endowments
        sim.e1[:] = np.random.exponential(par.zeta,size=par.N)
        sim.e2[:] = np.random.exponential(par.zeta,size=par.N)
        sim.e3[:] = np.random.exponential(par.zeta,size=par.N)

