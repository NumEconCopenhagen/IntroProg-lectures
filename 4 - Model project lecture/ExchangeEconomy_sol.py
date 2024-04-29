
# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import ipywidgets as widgets

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
        par.N = 50000
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
        sol = self.sol

        simvarnames = ['alphas','betas','e']
        for varname in simvarnames:
            sim.__dict__[varname] =  np.nan*np.ones((par.N,3)) 

    
        sim.x  =np.nan*np.ones((par.N,3))
        sim.I = np.nan*np.ones((par.N,1))


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
        sim.e[:] = np.random.exponential(par.zeta,size=(par.N,3))


    def plot_budget_shares(self):
        par = self.par 
        sim = self.sim

        
        fig,ax  = plt.subplots()
        for i in range(3):
            ax.hist(sim.betas[:,i],bins=100,alpha=0.7,density=True,label=f'Good {i+1}')


        ax.legend()

    def excess_demand(self,p):
        par = self.par
        sim = self.sim
        
        p_vec = np.array([p[0],p[1],1])

        sim.I[:,:] = (sim.e*p_vec).sum(axis=1).reshape((par.N,1))

        sim.x[:,:] = sim.betas*sim.I/p_vec

        demand = sim.x.sum(axis=0)

        supply = sim.e.sum(axis=0)


        return demand-supply
    
    def utility(self):
        par = self.par
        sim = self.sim

        u = np.prod(sim.x **sim.betas,axis=1)**par.gamma
        
        return u


    def find_equilibrium_scipy(self):
        par = self.par 
        sim = self.sim

        betas_mean = sim.betas.mean(axis=0)
        guess = np.array(betas_mean[:2]/betas_mean[2])
        print(f'Initial guess:\n P1 = {guess[0]:.2f}, P2 = {guess[1]:.2f}')
        
        obj = lambda p: self.excess_demand(p)[:2]
        sol = optimize.root(obj,guess)

        p = sol.x
        print(f'Solution:\n P1 = {p[0]:.2f}, P2 = {p[1]:.2f}')
        print(f'Excess demand: {self.excess_demand(p)}')

        return sol


    def find_equilibrium(self,maxiter=1000,printit=True):
        par = self.par 
        sim = self.sim
        sol = self.sol

        betas_mean = sim.betas.mean(axis=0)
        guess = np.array(betas_mean[:2]/betas_mean[2])
        if printit: print(f'Initial guess:\n P1 = {guess[0]:.2f}, P2 = {guess[1]:.2f}')

        p = guess
        for it in range(maxiter):
            excess = self.excess_demand(p)
            if np.abs(excess[:2]).max()<par.epsilon:
                if printit: print(f'Solution:\n P1 = {p[0]:.2f}, P2 = {p[1]:.2f}')
                if printit: print(f'Excess demand: {excess}')

                sol.p = p
                break
                
            if printit: print(f'Iteration {it+1}: p = [{p[0]:3.3f},{p[1]:3.2f}] -> Excess demand = {excess[:2]}')
            
            
            #p = p + par.kappa*excess[:2]/par.N


            # A faster updating rule taking into account, the excess demand for good 3
            p = p + par.kappa*(excess[:2]-excess[2]/2)/par.N


    def plot_utility_distribution(self,gamma=None):
        par = self.par 
        sim = self.sim

        if gamma is not None:
            gamma_org = par.gamma
            par.gamma = gamma       

        # Calculate utility in the equilibrium
        self.find_equilibrium(printit=False)
        utility = self.utility()

         # Calculate and print the mean and variance of utility
        mean= np.mean(utility)
        var = np.var(utility)


        # Plot the distribution of utility
        fig, ax = plt.subplots()
        ax.hist(utility, bins=30, edgecolor='black',label = fr'$\mu$= {mean:.3f}, $\sigma$= {var:.3f}')
        


        e_org = sim.e.copy()
        # Caclucalte if all endowments are equal
        sim.e[:,:] = sim.e.mean(axis=0)
        self.find_equilibrium(printit=False)
        utility = self.utility()
        mean= np.mean(utility)
        var = np.var(utility)
        ax.hist(utility, bins=30, edgecolor='black',label = fr'Equal e: $\mu$= {mean:.3f}, $\sigma$= {var:.3f}')
        
        
        ax.set_title('Distribution of utility')
        ax.set_xlabel('Utility')
        ax.set_ylabel('Number of agents')
        ax.legend()

        # reset
        if gamma is not None:
            par.gamma = gamma_org
        
        sim.e[:] = e_org



    def plot_utility_distribution_interactive(self):
        widgets.interact(self.plot_utility_distribution, gamma=widgets.FloatSlider(value=0.8, min=0.1, max=2.0, step=0.1))

    def plot_excess_demand(self):
        par = self.par
        sim = self.sim

        p1_values = np.linspace(0.5, 10, 100)
        p2_values = [0.5, 1.5, 2, 3]

        fig, ax = plt.subplots()

        colors_green = plt.cm.Greens(np.linspace(0.5, 1, len(p2_values)))
        colors_red = plt.cm.Reds(np.linspace(0.5, 1, len(p2_values)))

        for idx, p2 in enumerate(p2_values):
            excess_demand1 = []
            excess_demand2 = []
            for p1 in p1_values:
                excess_vec = self.excess_demand([p1,p2])
                excess_demand1.append(excess_vec[0])
                excess_demand2.append(excess_vec[1])

            ax.plot(p1_values, excess_demand1, color=colors_green[idx], label=f'P2={p2}, Good 1')
            ax.plot(p1_values, excess_demand2, color=colors_red[idx], label=f'P2={p2}, Good 2')
        
        ax.axhline(0, color='black', linewidth=0.5)
        ax.set_xlabel('P1')
        ax.set_ylabel('Excess Demand')
        ax.set_title('Excess Demand for Different Price Levels')
        ax.legend()
        plt.show()

    


    def _plot_excess_demand_interative(self, p2=0.5):
        par = self.par
        sim = self.sim

        p1_values = np.linspace(0.5, 10, 100)

        excess_demand1 = []
        excess_demand2 = []
        for p1 in p1_values:
            excess_vec = self.excess_demand([p1,p2])
            excess_demand1.append(excess_vec[0])
            excess_demand2.append(excess_vec[1])
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(p1_values, excess_demand1, color='green', label='Good 1')
        ax.plot(p1_values, excess_demand2, color='red', label='Good 2')
        ax.axhline(0, color='black', linewidth=0.5)
        ax.set_xlabel('P1')
        ax.set_ylabel('Excess Demand')
        ax.set_title('Excess Demand for Different Price Levels')
        ax.legend()

        plt.show()


    def plot_excess_demand_interactive(self):
        widgets.interact(self._plot_excess_demand_interative, p2=widgets.FloatSlider(value=1.0, min=0.5, max=3.0, step=0.1))


    def plot_excess_demand_3d(self):
        par = self.par
        sim = self.sim

        p1_values = np.linspace(0.5, 3, 100)
        p2_values = np.linspace(0.5, 3, 100)

        P1, P2 = np.meshgrid(p1_values, p2_values)

        excess_demand1 = np.empty_like(P1)
        excess_demand2 = np.empty_like(P1)

        for i in range(P1.shape[0]):
            for j in range(P1.shape[1]):
                p1 = P1[i, j]
                p2 = P2[i, j]

                excess_vec = self.excess_demand([p1,p2])
                excess_demand1[i, j] = excess_vec[0]
                excess_demand2[i, j] = excess_vec[1]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(P1, P2, excess_demand1, cmap='viridis')
        

        ax.set_xlabel('P1')
        ax.set_ylabel('P2')
        ax.set_zlabel('Excess Demand, good 1')
        ax.set_title('Excess Demand for Different Price Levels')
        plt.show()


        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(P1, P2, excess_demand2, cmap='viridis')
        
        ax.set_xlabel('P1')
        ax.set_ylabel('P2')
        ax.set_zlabel('Excess Demand, good 2')
        ax.set_title('Excess Demand for Different Price Levels')
        plt.show()


