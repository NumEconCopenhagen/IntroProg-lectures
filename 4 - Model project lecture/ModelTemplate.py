### This is a template/example of how to structure the .py-file for a model


# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})




# Define the model
class ModelName():

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
        par.alpha = 1.

        # Simulation options
        par.T = 100



    def allocate(self):
        '''
        Allocate arrays for simulation
        '''


        par = self.par
        sim = self.sim

        simvarnames = ['epsilon'] # Which arrays should be available in the simulation namespace

        for varname in simvarnames:
            sim.__dict__[varname] =  np.nan*np.ones(par.T) # Allocate the size of the arrays
        # In this function all arrays have the same size, this is not always the case


    