from types import SimpleNamespace
from copy import deepcopy

import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})

class WorkerClass:

    def __init__(self):
        """ initialize the worker class with default parameters """
        
        par = self.par = SimpleNamespace()

        par.alpha = 0.50
        par.kappa = 1.0
        par.nu = 1.0/(2*16.0**2)
        par.w = 1.0
        par.tau = 0.30

        par.varepsilon = 1.0
        par.sigma = 1.001
        par.rho = 1.001
        
        # set simple case as default
        self.utility = self.utility_cb
        self.L_opt = self.L_opt_cb

    def copy(self):

        other = WorkerClass()
        other.par = deepcopy(self.par)
        other.utility = self.utility
        other.L_opt = self.L_opt

        return other
    
    def utility_cb(self,L,G):
        """ Cobb-Douglas utility function"""

        par = self.par

        C = par.kappa+(1-par.tau)*par.w*L
        CG = C**par.alpha*G**(1-par.alpha)

        return np.log(CG) - par.nu*L**(1+par.varepsilon)/(1+par.varepsilon)

    def L_opt_cb(self):
        """ solve the optimization problem analytically for Cobb-Douglas utility function """
        
        par = self.par

        wt = (1-par.tau)*par.w
        L_opt = (-par.kappa + np.sqrt(par.kappa**2+4*par.alpha/par.nu*wt**2))/(2*wt)

        return L_opt

    def utility_ces(self,L,G):
        """ CES utility function """

        par = self.par
        
        C = par.kappa+(1-par.tau)*par.w*L
        sigma_par = (par.sigma-1)/par.sigma
        ces = ( par.alpha* C**(sigma_par) + (1-par.alpha) *G**(sigma_par) )**(1/sigma_par)
        
        return  ( (ces)**(1-par.rho) -1) /(1-par.rho) - par.nu*L**(1+par.varepsilon)/(1+par.varepsilon)
    
    def L_opt_ces(self,G=1):
        """ solve the optimization problem numerically for CES utility function """
        
        return self.L_opt_general('utility_ces',G)
    
    def L_opt_general(self,utility_func_name,G=1.0):
        """ solve the optimization problem numerically for a given utility function """

        par = self.par

        def obj(L):
            return - getattr(self,utility_func_name)(L,G)

        L_opt = optimize.minimize_scalar(obj,bounds=(1e-8,24-1e-8),method='bounded').x

        return L_opt
    
    def test_L_opt(self,G=1.0, utility_func_name='utility_cb'):
        """ test the L_opt function """

        L_opt = self.L_opt_general(utility_func_name,G=G)
        L_opt_analytical = self.L_opt_cb()
        L_grid = np.linspace(np.fmax(1e-8,L_opt-1),np.fmin(24-1e-8,L_opt+1),100)
        us = np.empty(L_grid.size)

        for i,L in enumerate(L_grid):
            us[i] = getattr(self,utility_func_name)(L,G)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title(f'G={G:.2f}')
        ax.plot(L_grid,us,label='Utility')
        ax.axvline(L_opt_analytical,color='green',linestyle='solid',label='Analytical solution')
        ax.axvline(L_opt,color='red',linestyle='dashed', label = 'Numerical solution')
        ax.legend()
        ax.set_xlabel('L')
        ax.set_ylabel('utility')

    def plot_L_opt_w(self,G=1.0,L_opt_name='L_opt'):
        """Plot optimal L across wages"""
        par = self.par

        w_grid = np.linspace(0.1,10.0,100)
        L_opts = np.empty(w_grid.size)

        original_w = par.w
        for i,w in enumerate(w_grid):
            par.w = w
            L_opts[i] =  getattr(self,L_opt_name)()

        par.w = original_w # restore

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.set_title(f'G={G:.2f}')
        ax.plot(w_grid,L_opts)
        ax.set_xlabel('w')
        ax.set_ylabel('optimal L')

class GovClass:

    def __init__(self):
        """ initialize the government class with default parameters"""

        self.sol_cb = SimpleNamespace() # Solutions for Cobb-Douglas worker
        self.sol_ces = SimpleNamespace() # Solutions for CES worker
        self.worker = WorkerClass()
        
    def copy(self):

        other = GovClass()
        other.sol_cb = deepcopy(self.sol_cb)
        other.sol_ces = deepcopy(self.sol_ces)
        other.worker = self.worker.copy()

        return other

    def social_value(self,tau):
        """ compute the social value for a given tau assuming the worker is a Cobb-Douglas worker """

        worker = self.worker
        par = worker.par
        
        par.tau = tau
        L_opt = worker.L_opt()
        G_opt = par.tau*par.w*L_opt

        return worker.utility(L_opt,G_opt),L_opt,G_opt

    def plot_social_value_tau(self):
        """ plot the social value as a function of tau assuming the worker is a Cobb-Douglas worker"""

        tau_grid = np.linspace(0.01,0.99,100)
        sv = np.empty(tau_grid.size)
        L_opt = np.empty(tau_grid.size)
        G_opt = np.empty(tau_grid.size)

        for i,tau in enumerate(tau_grid):
            sv[i],L_opt[i],G_opt[i] = self.social_value(tau)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(tau_grid,sv)
        ax.set_xlabel('$\\tau$')
        ax.set_ylabel('social value')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(tau_grid,L_opt)
        ax.set_xlabel('$\\tau$')
        ax.set_ylabel('L')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(tau_grid,G_opt)
        ax.set_xlabel('$\\tau$')
        ax.set_ylabel('G')

    def find_optimal_tau(self,type='cb'):
        """ Find the optimal tau for a given worker type"""

        sol = getattr(self,'sol_'+type)
        
        worker = self.worker
        par = worker.par

        if type == 'cb':

            def objective(tau):
                par.tau = tau
                L_opt = worker.L_opt_cb()
                G_opt = par.tau*par.w*L_opt
                return -worker.utility_cb(L_opt,G_opt)

        elif type == 'ces':
            
            def objective(tau):
                par.tau = tau
                self.find_G(tau,plotit=False)
                return -worker.utility_ces(sol.L_opt,sol.G_opt)
            
        res = optimize.minimize_scalar(objective,bounds=(0.01,0.99),method='bounded')
        sol.tau_opt  = res.x

        if type == 'cb':
            sol.sv_opt,sol.L_opt,sol.G_opt = self.social_value(sol.tau_opt)
        elif type == 'ces':
            self.find_G(sol.tau_opt,plotit=False)
            sol.sv_opt = -res.fun
            

        print(f'tau_opt = {sol.tau_opt:.2f}')
        print(f'social value = {sol.sv_opt:.2f}')
        print(f'L_opt = {sol.L_opt:.2f}')
        print(f'G_opt = {sol.G_opt:.2f}')

        # illustrate optimal tau
        tau_grid = np.linspace(0.01,0.99,100)
        sv = np.empty(tau_grid.size)

        for i,tau in enumerate(tau_grid):
            if type == 'cb':
                sv[i],_,_ = self.social_value(tau)
            elif type == 'ces':
                self.find_G(tau,plotit=False)
                sv[i] = worker.utility_ces(sol.L_opt,sol.G_opt) 
            
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(tau_grid,sv)
        ax.axvline(sol.tau_opt,color='red',linestyle='dashed')
        ax.set_xlabel('$\\tau$')       
        ax.set_ylabel('sociale value')       

    def plot_budget_constriant(self,tau):
        """ plot budget constraint for given tau, assuming worker is CES"""

        sol_ces = self.sol_ces

        worker = self.worker
        par = worker.par

        par.tau = tau

        G_grid = np.linspace(0.01,10.0,100)
        BC_grid = np.empty(G_grid.size)
        L_opt_grid = np.empty(G_grid.size)

        for i,G in enumerate(G_grid):
            L_opt_grid[i] = worker.L_opt_ces(G)
            BC_grid[i] = par.tau*par.w*L_opt_grid[i] - G
        
        fig = plt.figure(figsize=(10,4))

        ax1 = fig.add_subplot(1,2,1)
        ax1.plot(G_grid,L_opt_grid)
        ax1.axvline(sol_ces.G_opt,color='red',linestyle='dashed')
        ax1.set_xlabel('G')
        ax1.set_ylabel('labor supply')
        
        ax2 = fig.add_subplot(1,2,2)
        ax2.plot(G_grid,BC_grid)
        ax2.axvline(sol_ces.G_opt,color='red',linestyle='dashed')
        ax2.set_xlabel('G')
        ax2.set_ylabel('budget surplus')
        fig.tight_layout()
        
    def find_G(self,tau,plotit=True):
        """ find G for CES formulation """

        sol_ces = self.sol_ces
        sol_cb = self.sol_cb

        worker = self.worker
        par = worker.par

        def obj(G):
            par.tau = tau 
            L_opt = worker.L_opt_ces(G)
            
            return par.tau*par.w*L_opt - G
        
        opt = optimize.root_scalar(obj,x0 =sol_cb.G_opt ,x1 = sol_cb.G_opt*0.9 )
        
        sol_ces.G_opt = opt.root
        sol_ces.L_opt = worker.L_opt_ces(sol_ces.G_opt)

        if plotit: self.plot_budget_constriant(tau)