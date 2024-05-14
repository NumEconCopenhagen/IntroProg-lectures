
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



    def u_cb(self,L, G):
        '''
        Cobb Douglass utility function
        '''
        par = self.par
        C = par.kappa + (1-par.tau)*par.w*L
        C_bounded = np.fmax(C,par.zero_tol)
        G_bounded = np.fmax(G,par.zero_tol)

        return np.log( (C_bounded**par.alpha) * G_bounded**(1-par.alpha) ) - par.nu*L**2/2
    

    def L_star_analytical(self,G=1):
        '''
        Analytical solution for optimal labor supply
        '''
        par = self.par
        sol_cb = self.sol_cb

        wtilde = par.w*(1-par.tau)

        sol_cb.L_analytical = (-par.kappa+ np.sqrt(par.kappa**2+ 4*par.alpha/par.nu*wtilde**2)) /(2*wtilde)
     
    def u_ces(self,L,G):
        '''
        Utility function
        '''
        par = self.par
        C = par.kappa + (1-par.tau)*par.w*L
        C_bounded = np.fmax(C,par.zero_tol)
        G_bounded = np.fmax(G,par.zero_tol)

        term1 = (par.alpha * C_bounded**((par.sigma-1)/par.sigma) + (1-par.alpha) * G_bounded**((par.sigma-1)/par.sigma))**(par.sigma/(par.sigma-1))
        term2 = par.nu * L**(1+par.varepsilon) / (1+par.varepsilon)

        return (term1**(1-par.rho) -1) / (1-par.rho) - term2
    