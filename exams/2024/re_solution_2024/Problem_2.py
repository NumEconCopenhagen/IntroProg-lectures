import numpy as np
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scipy import optimize

class FSModel:
    
    def __init__(self,**kwargs):
    
        '''
        Initialize the model with parameters
        '''
        
        par = self.par = SimpleNamespace()
        self.sim = SimpleNamespace()

        # deterministic
        par.Y0 = 1
        par.B0 = 0
        par.g = 0.02
        par.r = 0.04
        par.PD = 0.1
        par.T = 500

        # stochastic
        par.mug = 0.02
        par.sigmag = 0.02
        par.mur = 0.00
        par.sigmar = 0.02
        par.nu = 0.001
        par.K = 1000

        # extended
        par.xi = 0.1
        par.kappa1 = 0.01 
        par.kappa2 = 0.5

        for k,v in kwargs.items():
            setattr(par,k,v)

        self.allocate()

    def allocate(self):
        '''
        Allocate size of simulation arrays
        '''

        par = self.par
        sim = self.sim

        sim.Y = np.zeros((par.T,par.K))
        sim.B = np.zeros((par.T,par.K))
        sim.b = np.zeros((par.T,par.K))

        sim.g = np.zeros((par.T,par.K))
        sim.r = np.zeros((par.T,par.K))

        sim.epsg = np.zeros((par.T,par.K))
        sim.epsr = np.zeros((par.T,par.K))

        sim.check = np.zeros(par.K,dtype=bool)
       
    def simulate(self):
        '''
        Simulate the non-stochastic model
        
        '''

        par = self.par
        sim = self.sim

        sim.Y[0,0] = par.Y0
        sim.B[0,0] = par.B0
        sim.b[0,0] = par.B0/par.Y0

        for t in range(1,par.T):

            sim.Y[t,0] = sim.Y[t-1,0]*(1+par.g)
            sim.B[t,0] = sim.B[t-1,0]*(1+par.r) + par.PD*sim.Y[t,0]
            sim.b[t,0] = sim.B[t,0]/sim.Y[t,0]

    def simulate_stochastic(self,rng=None,max_value=2.0,extended=False):
        '''
        Simulate the stochastic model
        '''

        par = self.par
        sim = self.sim
        
        if not rng is None:
            sim.epsg[:] =  rng.normal(par.mug,par.sigmag,(par.T,par.K)) 
            sim.epsr[:] =  rng.normal(par.mur,par.sigmar,(par.T,par.K))

        sim.Y[0,:] = par.Y0
        sim.B[0,:] = par.B0
        sim.b[0,:] = par.B0/par.Y0
        sim.r[0,:] = np.exp(sim.epsr[0,:])-1+par.nu*sim.b[0,:]**2

        if extended:
            sim.g[0,:] = np.exp( sim.epsg[0,:])-1-par.xi*(sim.r[0,:]-par.mur ) + par.kappa1 * par.PD**par.kappa2  
        else:
            sim.g[0,:] = np.exp(sim.epsg[0,:])-1


        sim.check[:] = True        
        for t in range(1,par.T):

            I = sim.check

            sim.Y[t,I] = (1+sim.g[t-1,I])*sim.Y[t-1,I]
            sim.B[t,I] = (1+ sim.r[t-1,I])*sim.B[t-1,I] + par.PD*sim.Y[t,I]
            sim.r[t,I] = np.exp( sim.epsr[t,I])-1+par.nu*sim.b[t,I]**2 
            
            if extended:
                sim.g[t,I] = np.exp( sim.epsg[t,I])-1-par.xi*(sim.r[t,I]-par.mur ) + par.kappa1*par.PD**par.kappa2  
            else:
                sim.g[t,I] = np.exp( sim.epsg[t,I])-1 

            sim.b[t,I] = sim.B[t,I]/sim.Y[t,I]

            sim.check[I] = sim.check[I] & (sim.b[t,I] < max_value)

            if np.all(sim.check==False): break

    def find_PD_limit_max_value(self,PD_interval=[0,1],max_value=2.0):
        '''
        Find the maximum PD that still allows for debt at some value max_value
        '''

        par = self.par
        sim = self.sim

        def obj(PD):
            
            par.PD = PD
            self.simulate()
            return np.max(sim.b[:,0]) - max_value
        
        res = optimize.root_scalar(obj,bracket=PD_interval)
        return res.root
    
    def plot_b(self,i=0):
        '''
        Plot the debt ratio over time 
        '''

        par = self.par
        sim = self.sim
        fig, ax = plt.subplots(1,1)

        ax.plot(range(par.T),sim.b[:,i])
        ax.set_xlabel('T')
        ax.set_ylabel('B/Y')
        plt.show()

    def plot_sustainable_probability(self,PD_interval=[0.0,0.10],N=100,
                                     extended=False,
                                     ax=None,show= True,label=None):
        '''
        Plot the probability of sustainable debt as a function of PD    
        '''

        par = self.par
        sim = self.sim

        if ax is None:
            fig, ax = plt.subplots(1,1)

        PD_grid = np.linspace(PD_interval[0],PD_interval[1],N)
        s_prob = np.zeros(N)

        for i,PD in enumerate(PD_grid):
            par.PD = PD
            self.simulate_stochastic(extended=extended)
            s_prob[i] = np.mean(sim.check)

        ax.plot(PD_grid, s_prob,label=label)
        ax.set_xlabel('PD')
        ax.set_ylabel('propability of sustainable debt')
        if show: plt.show()

        return ax