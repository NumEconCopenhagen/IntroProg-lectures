
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})




# Define the model
class Hairdresser():

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

        # Draw shocks
        self.simulate()


    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par

        # Model parameters
        par.eta = 0.5
        par.w = 1.0
        par.rho = 0.9

        par.iota = 0.01
        par.sigma_epsilon = 0.1

        par.R =(1+0.01)**(1/12)


        # Simulation options
        par.T = 120
        par.K = 10000



    def allocate(self):
        '''
        Allocate arrays for simulation
        '''


        par = self.par
        sim = self.sim

        simvarnames = ['epsilon','log_kappa','kappa','ell','profits'] # Which arrays should be available in the simulation namespace

        for varname in simvarnames:
            sim.__dict__[varname] =  np.nan*np.ones((par.K,par.T)) # Allocate the size of the arrays
        
        sim.h = np.nan*np.ones((par.K)) # Allocate the size of the arrays


    
    def simulate(self,seed=1234):
        '''
        Draw random shocks for your model

        It's a good idea to draw all the random shocks you want in one function and store them 
        You can then be explicit about when you change the seed and take a new draw.
        '''
        par = self.par
        sim = self.sim 
        
        np.random.seed(seed)

        sim.epsilon[:] = np.random.normal(loc= -0.5*par.sigma_epsilon**2,scale=par.sigma_epsilon, size=(par.K,par.T))
        
        
        sim.log_kappa[:,0] = 0+ sim.epsilon[:,0]
        for t in range(1,par.T):
            sim.log_kappa[:,t] = par.rho*sim.log_kappa[:,t-1] + sim.epsilon[:,t]

        sim.kappa[:,:] = np.exp(sim.log_kappa)

    def H_func(self, policy='static'):
        par = self.par
        sim = self.sim
        sol = self.sol

        # calculate h
        if policy == 'static':
            sim.ell[:,:] = self.ell_star_analytical(sim.kappa)
        else:
            assert False, 'Not implemented yet'


        sim.profits[:,0] = self.profits(sim.kappa[:,0],sim.ell[:,0],0)
        for t in range(1,par.T):
            sim.profits[:,t] = self.profits(sim.kappa[:,t],sim.ell[:,t],sim.ell[:,t-1])
        
        sim.h[:] = np.sum(  par.R**(-np.arange(0,par.T,1))* sim.profits[:,:] ,axis=1)

        sol.H = np.mean(sim.h)


        

    def profits(self,kappa,ell,ell_lag):
        '''
        Compute profits
        '''
        par = self.par

        y = np.fmax(ell,1e-8)
        
        return self.profits_static(kappa,y) - (~np.isclose(ell,ell_lag))* par.iota

    def profits_static(self,kappa,ell):
        '''
        Compute _static
        '''
        par = self.par

        return kappa*ell**(1-par.eta) - par.w*ell

    def ell_star_analytical(self,kappa):
        '''
        Compute analytical solution for optimal labor
        '''
        par = self.par

        return ((1-par.eta)*kappa/par.w)**(1/par.eta)
    

    def plot_profits_across_ell(self,kappa):
        '''
        Plot profits as a function of labor
        '''
        par = self.par

        ell_grid = np.linspace(0.1,5,100)
        profits = self.profits_static(kappa,ell_grid)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(ell_grid,profits,label=r'$\Pi(\ell)$')
        ax.set_xlabel('$\ell$')
        ax.set_ylabel('profits')

        ax.axvline(self.ell_star_analytical(kappa),color='red',linestyle='-',label='$\ell^*$')

        # Add numerical solution

        ell_star = optimize.minimize_scalar(lambda ell: -self.profits_static(kappa,ell)).x
        ax.axvline(ell_star,color='green',linestyle='--',label='$\ell^*_{num}$')

        ax.legend()