
# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})




# Define the model
class GovModel():

    def __init__(self, **kwargs):
        '''
        Initialize the model with default parameters
        kwargs allow any parameter in the par namespace to be overridden by the user
        '''

        self.par = par = SimpleNamespace() # Create a namespace object for parameters
        self.sol_cb = sol_cb = SimpleNamespace() # Create a namespace object for solution results
        self.sol_ces = sol_ces = SimpleNamespace() # Create a namespace object for solution results

        # Set default parameters
        self.setup()

        # Update parameters with user input
        for key, value in kwargs.items():
            setattr(par, key, value)


    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par

        # Model parameters
        par.alpha = 0.5
        par.kappa =1.0
        par.nu = 1/(2*16**2)
        par.w = 1.0
        par.tau = 0.3

        par.sigma = 1.001 
        par.rho = 1.001
        par.varepsilon = 1.0
        
        par.zero_tol = 1e-8