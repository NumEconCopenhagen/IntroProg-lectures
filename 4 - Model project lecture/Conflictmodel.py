
# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize, interpolate
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})


# Define the model
class ConflictModel():

    def __init__(self, **kwargs):
        '''
        Initialize the model with default parameters
        kwargs allow any parameter in the par namespace to be overridden by the user
        '''

        self.par = par = SimpleNamespace() # Create a namespace object for parameters
        self.sol = sol = SimpleNamespace() # Create a namespace object for solution results
        

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
        par.ea = 10
        par.eb = 10
        par.epsilon = 10.
        par.eta = 1000. 

        par.small = 1e-4
