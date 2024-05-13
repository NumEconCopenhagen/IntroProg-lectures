
from types import SimpleNamespace
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})

class SalonClass:

    def __init__(self):
        """ initialize the Salon class with default parameters """

        par = self.par = SimpleNamespace()
        sim = self.sim = SimpleNamespace()

        # parameters for static version
        par.w = 1.0
        par.eta = 0.5
        
        # parameters for dynamic version
        par.rho = 0.90
        par.R = (1.0+0.01)**(1/12)
        par.iota = 0.01
        par.sigma_epsilon = 0.10
        par.K = 10_000
        par.T = 120

        par.Delta = 0.05

        par.N = 10 # number of retakes for ensuring optimal Delta is found
        par.adjust = 0.98 # adjustment factor optimal ell 
        par.Delta_adjust = -0.05 # adjustment for delta across t

    def static_ell_opt(self,kappa):
        """ optimal labor demand in static model """

        par = self.par
        return ((1-par.eta)*kappa/par.w)**(1/par.eta)

    def p(self,kappa,y):
        """ price as a function of demand """

        par = self.par
        return kappa*y**(-par.eta)
                
    def profits_static(self,kappa,ell):
        """ static profits"""

        par = self.par

        y = ell
        p = self.p(kappa,y)

        return p*y - par.w*ell

    def plot_profits_static(self,kappa):
        """ plot profits in the static model"""
        
        par = self.par

        ell_grid = np.linspace(0.01,2.0,100)
        profits = self.profits_static(kappa,ell_grid)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title(f'$\\kappa={kappa:.2f}$')
        ax.plot(ell_grid,profits,label='profits')
        ell_opt = self.static_ell_opt(kappa)
        ax.axvline(ell_opt,color='red',label='$\ell^*$')
        ax.set_xlabel('ell')
        ax.set_ylabel('profits')


        # Add numerical solution

        ell_star = optimize.minimize_scalar(lambda ell: -self.profits_static(kappa,ell), bounds =(0.01,2)).x
        ax.axvline(ell_star,color='green',linestyle='--',label='$\ell^*_{num}$')
        ax.legend()
        


    def simulate(self,state=None):
        """ simulate the dynamic model"""
        
        if not state is None:
            np.random.set_state(state)

        sim = self.sim
        par = self.par

        # preallocate simulation 
        sim.epsilon = np.empty((par.K,par.T))
        sim.kappa = np.empty((par.K,par.T))
        sim.profits = np.empty((par.K,par.T))
        sim.ell = np.empty((par.K,par.T))
        sim.h = np.empty(par.K)

        # Draw shocks 
        sim.epsilon[:,:] =  np.random.normal(-0.5*par.sigma_epsilon**2,par.sigma_epsilon,size=(par.K,par.T))        
        log_kappa = np.zeros((par.K,par.T))


        # Calculate log(kappa)
        log_kappa[:,0] = sim.epsilon[:,0]
        for t in range(1,par.T):
            log_kappa[:,t] = par.rho*log_kappa[:,t-1] + sim.epsilon[:,t] 
        
        # Calculate kappa
        sim.kappa[:,:] = np.exp(log_kappa)
        
    def profits(self,kappa,ell,ell_lag):
        """ profits in the dynamic model"""
    
        par = self.par

        y = np.fmax(ell,1e-8)
        p = self.p(kappa,y)

        return p*y - par.w*ell - (~np.isclose(ell,ell_lag))*par.iota

    def evaluate(self,policy='static_policy',vectorized=True,print_result=True):
        """ evaluate a given policy in the dynamic model """
        
        par = self.par
        sim = self.sim

        rule = getattr(self,policy)
        
        if not vectorized:

            # not utilizing vectorization across k
            for k in range(par.K):

                ell_lag = 0 # initial 
                
                for t in range(par.T):

                    kappa = sim.kappa[k,t]
                    ell = rule(kappa,args={'ell_lag':ell_lag,'t':t})
                    sim.ell[k,t] = ell
                    sim.profits[k,t] = self.profits(kappa,ell,ell_lag)
                    ell_lag = ell
        
                sim.h[k] = np.sum( par.R**(-np.arange(0,par.T,1))*sim.profits[k,:])

        elif vectorized:

            ell_lag = 0 # initial 
            for t in range(par.T):
                kappa = sim.kappa[:,t]
                ell = rule(kappa,args={'ell_lag':ell_lag,'t':t})
                sim.ell[:,t] = ell
                sim.profits[:,t] = self.profits(kappa,ell,ell_lag)
                ell_lag = ell
            
            sim.h[:] = np.sum(par.R**(-np.arange(0,par.T,1))*sim.profits[:,:],axis=1)

        sim.H = np.mean(sim.h)

        if print_result:
            print(f'H = {sim.H:.2f} for {policy}')
            if 'delta' in policy:
                print(f'Delta = {par.Delta:.4f}')

    def static_policy(self,kappa,args={}):
        """ static policy function"""

        return self.static_ell_opt(kappa)
    
    def delta_policy(self,kappa,args={}):
        """ delta policy function"""

        ell_lag = args['ell_lag']
        ell_star = self.static_ell_opt(kappa)
        shift = (np.abs(ell_star - ell_lag)>self.par.Delta)        
        return shift*ell_star + (1-shift)*ell_lag
    
    def alt_delta_policy(self,kappa,args={}):
        """ 
        Alternative Delta policy function where optimal ell is 
        allowed to be different from analytical solution
        """
        
        par = self.par
        
        ell_lag = args['ell_lag']
        ell_star = self.static_ell_opt(kappa)
        shift = ( np.abs(ell_star*par.adjust - ell_lag)  > par.Delta)    

        return shift*ell_star*par.adjust + (1-shift)*ell_lag

    def alt2_delta_policy(self,kappa,args={}):
        """ 
        Alternative Delta policy function where optimal ell is allowed to be different 
        from analytical solution and the adjustment is time dependent 
        """

        par = self.par
        ell_lag = args['ell_lag']
        t = args['t']
        ell_star = self.static_ell_opt(kappa)
        shift = ( np.abs(ell_star*par.adjust - ell_lag) > ((par.T-t)/par.T)**(par.Delta_adjust ) * par.Delta)    

        return shift*ell_star*par.adjust + (1-shift)*ell_lag

    def find_optimal_pars(self,policy='delta_policy',opt_pars=['Delta'],bounds=((0,None),),print_result=True):
        """ find optimal parameters for a given policy and simulation """
        
        par = self.par
        x0 = [getattr(par,par_name) for par_name in opt_pars]

        def obj(x):
        
            for i,par_name in enumerate(opt_pars):
                setattr(par,par_name,x[i])
            
            self.evaluate(policy=policy,print_result=False)
            return -self.sim.H

        res = optimize.minimize(obj,x0=x0,bounds=bounds,method='Nelder-Mead',tol=1e-8)
        
        if print_result:
            for par_name in opt_pars:
                print(f'Optimized {par_name} = {getattr(par,par_name):.4f}')
            print(f'Optimal H       = {self.sim.H:.4f}')

    def ensure_optimal_pars(self,policy ='delta_policy',opt_pars = ['Delta'],bounds=((0,None),),state=None):
        """ Ensure that the optimal parameters are found for a given policy by 
        repeating process for par.N different simulations"""
        
        if not state is None:
            np.random.set_state(state)

        par = self.par
        sim = self.sim

        sim.optimal_delta_grid = np.empty(par.N)
        for par_name in opt_pars:
            setattr(sim,f'optimal_{par_name}_grid',np.empty(par.N))

        sim.optimal_H_grid = np.empty(par.N)

        for i in range(par.N):

            self.simulate(state=None)
            self.find_optimal_pars(policy=policy,opt_pars=opt_pars, bounds =bounds,print_result=False)
            
            for par_name in opt_pars:
                getattr(sim,f'optimal_{par_name}_grid')[i] = getattr(par,par_name)

            sim.optimal_H_grid[i] = sim.H

        for par_name in opt_pars:
            par_grid = getattr(sim,f'optimal_{par_name}_grid')
            print(f'E({par_name})   = {np.mean(par_grid):.4f}')
            print(f'Std({par_name}) = {np.std(par_grid):.4f}')

        print(f'E(H)       = {np.mean(sim.optimal_H_grid):.4f}')

    def avg_adjustments(self):
        """ print average adjustments in given solution"""

        print(f'{np.isclose(self.sim.ell[:,1:],self.sim.ell[:,:-1]).mean():.4f}')

#### Plotting functions ####

    def plot_kappa(self,k=0):

        sim = self.sim
        fig = plt.figure(figsize=(5,4))
        ax = fig.add_subplot(1,1,1)
        ax.plot(list(range(self.par.T)), sim.kappa[k,:])
        ax.set_xlabel('t')
        ax.set_ylabel(r'$\kappa_t$')

    def plot_kappa_agg(self,k=0):
        sim = self.sim
        fig = plt.figure(figsize=(5,4))
        ax = fig.add_subplot(1,1,1)
        ax.plot(list(range(self.par.T)), sim.kappa.mean(axis=0),label='mean')
        ax.plot(list(range(self.par.T)), np.median(sim.kappa,axis=0),label='median')
        ax.set_xlabel('t')
        ax.set_ylabel(r'$ \kappa_t $')
        ax.legend()

    def plot_policy(self,k=0,plot_kappa=True):

        sim = self.sim

        fig = plt.figure(figsize=(10,4))
        ax1 = fig.add_subplot(1,2,1)
        ax1.plot(list(range(self.par.T)), sim.ell[k,:])
        ax1.set_xlabel('T')
        ax1.set_ylabel(r'$\ell$')

        ax2 = fig.add_subplot(1,2,2)
        ax2.plot(list(range(self.par.T)), sim.profits[k,:])
        ax2.set_xlabel('t')
        ax2.set_ylabel('profits')

        fig.tight_layout()
        
        if plot_kappa:
            self.plot_kappa(k=k)

    def plot_policy_agg(self,plot_kappa=True):

        sim = self.sim

        fig = plt.figure(figsize=(10,4))

        ax1 = fig.add_subplot(1,2,1)
        ax1.plot(list(range(self.par.T)), sim.ell.mean(axis=0),label='mean')
        ax1.plot(list(range(self.par.T)), np.median(sim.ell,axis=0),label='median')
        ax1.set_xlabel('T')
        ax1.set_ylabel(r'$\ell$')
        ax1.legend()
        
        ax2 = fig.add_subplot(1,2,2)
        ax2.plot(list(range(self.par.T)), sim.profits.mean(axis=0),label='mean')
        ax2.plot(list(range(self.par.T)), np.median(sim.profits,axis=0),label='median')
        ax2.set_xlabel('t')
        ax2.set_ylabel('profits')
        ax2.legend()

        fig.tight_layout()

        if plot_kappa:
            self.plot_kappa_agg()

    def plot_H_across_delta(self,policy='delta_policy'):

        opt_delta = self.par.Delta
        delta_grid = np.linspace(0.0,3.0,50)*opt_delta
        H_grid = np.empty(delta_grid.shape)

        for i, delta in enumerate(delta_grid):
            self.par.Delta = delta
            self.evaluate(policy=policy,print_result=False)
            H_grid[i] = self.sim.H
        
        # reset delta
        self.par.Delta = opt_delta

        fig = plt.figure(figsize=(5,4))
        
        ax = fig.add_subplot(1,1,1)
        ax.plot(delta_grid,H_grid)
        ax.axvline(opt_delta,color='red',linestyle='dashed')
        ax.set_xlabel(r'$\Delta$')
        ax.set_ylabel('H')
