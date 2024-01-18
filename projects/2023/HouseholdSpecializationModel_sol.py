
from types import SimpleNamespace

import numpy as np
from scipy import optimize

import pandas as pd 
import matplotlib.pyplot as plt

class HouseholdSpecializationModelClass:

    def __init__(self):
        """ setup model """

        # a. create namespaces
        par = self.par = SimpleNamespace()
        sol = self.sol = SimpleNamespace()

        # b. preferences
        par.rho = 2.0
        par.nu = 0.001
        par.epsilon = 1.0
        par.omega = 0.5 

        # c. household production
        par.alpha = 0.5
        par.sigma = 1.0

        # d. wages
        par.wM = 1.0
        par.wF = 1.0
        par.wF_vec = np.linspace(0.8,1.2,5)

        # e. targets
        par.beta0_target = 0.4
        par.beta1_target = -0.1

        # f. solution
        sol.LM_vec = np.zeros(par.wF_vec.size)
        sol.HM_vec = np.zeros(par.wF_vec.size)
        sol.LF_vec = np.zeros(par.wF_vec.size)
        sol.HF_vec = np.zeros(par.wF_vec.size)

        sol.beta0 = np.nan
        sol.beta1 = np.nan

    def calc_utility(self,LM,HM,LF,HF):
        """ calculate utility """

        par = self.par
        sol = self.sol

        # a. consumption of market goods
        C = par.wM*LM + par.wF*LF

        # b. home production
        
        # DELETE
        if np.isclose(par.sigma,0.0):    
            H = np.fmin(HM,HF)
        elif np.isclose(par.sigma,1.0):
            H = HM**(1-par.alpha)*HF**par.alpha
        else:
            _M_term = (1-par.alpha)*np.fmax(HM,1e-8)**((par.sigma-1)/par.sigma)
            _F_term = par.alpha*np.fmax(HF,1e-8)**((par.sigma-1)/par.sigma)
            _inner = _M_term + _F_term
            H = np.fmax(_inner,1e-8)**(par.sigma/(par.sigma-1))
        # REPLACE
        # H = HM**(1-par.alpha)*HF**par.alpha
        # END

        # c. total consumption utility
        Q = C**par.omega*H**(1-par.omega)
        utility = np.fmax(Q,1e-8)**(1-par.rho)/(1-par.rho)

        # d. disutlity of work
        epsilon_ = 1+1/par.epsilon
        TM = LM+HM
        TF = LF+HF
        disutility = par.nu*(TM**epsilon_/epsilon_+TF**epsilon_/epsilon_)
        
        return utility - disutility

    def solve_discrete(self,do_print=False):
        """ solve model discretely """
        
        par = self.par
        sol = self.sol
        opt = SimpleNamespace()
        
        # a. all possible choices
        x = np.linspace(0,24,49)
        LM,HM,LF,HF = np.meshgrid(x,x,x,x) # all combinations
    
        LM = LM.ravel() # vector
        HM = HM.ravel()
        LF = LF.ravel()
        HF = HF.ravel()

        # b. calculate utility
        u = self.calc_utility(LM,HM,LF,HF)
    
        # c. set to minus infinity if constraint is broken
        I = (LM+HM > 24) | (LF+HF > 24) # | is "or"
        u[I] = -np.inf
    
        # d. find maximizing argument
        j = np.argmax(u)
        
        opt.LM = LM[j]
        opt.HM = HM[j]
        opt.LF = LF[j]
        opt.HF = HF[j]

        # e. print
        if do_print:
            for k,v in opt.__dict__.items():
                print(f'{k} = {v:6.4f}')

        return opt

    def solve(self,do_print=False):
        """ solve model continously """

        # DELETE

        sol = self.sol
        opt = SimpleNamespace()

        # a. initial guess from discrete solution
        opt = self.solve_discrete()
        x0 = np.array([opt.LM,opt.HM,opt.LF,opt.HF])

        # b. run optimizer
        def obj(x):

            LM = x[0]
            HM = x[1]
            LF = x[2]
            HF = x[3]

            TM = LM+HM
            TF = LF+HF

            if TM > 24:
                LM *= 24/TM
                HM *= 24/TM
            
            if TF > 24:
                LF *= 24/TF
                HF *= 24/TF

            return -self.calc_utility(LM,HM,LF,HF)

        result = optimize.minimize(obj,x0,method='Nelder-Mead',bounds=((0,24),(0,24),(0,24),(0,24)))

        opt.LM = result.x[0]
        opt.HM = result.x[1]
        opt.LF = result.x[2]
        opt.HF = result.x[3]   

        # c. print
        if do_print:
            for k,v in opt.__dict__.items():
                print(f'{k} = {v:6.4f}')

        return opt
    
    def solve_wF_vec(self,discrete=False):
        """ solve model for vector of female wages """

        # DELTE

        par = self.par
        sol = self.sol

        for i,wF in enumerate(par.wF_vec):
        
            self.par.wF = wF

            if discrete:
                opt = self.solve_discrete()
            else:
                opt =  self.solve()

            sol.LM_vec[i] = opt.LM
            sol.HM_vec[i] = opt.HM
            sol.LF_vec[i] = opt.LF
            sol.HF_vec[i] = opt.HF

    def run_regression(self):
        """ run regression """

        par = self.par
        sol = self.sol

        x = np.log(par.wF_vec)
        y = np.log(sol.HF_vec/sol.HM_vec)
        A = np.vstack([np.ones(x.size),x]).T
        sol.beta0,sol.beta1 = np.linalg.lstsq(A,y,rcond=None)[0]
    
    def estimate(self,alpha=None,sigma=None):
        """ estimate alpha and sigma """

        # DELTE

        par = self.par
        sol = self.sol

        # a. objective
        def obj(x):

            par.alpha = x[0]
            par.sigma = x[1]

            # a. solve
            self.solve_wF_vec()

            # b. run regression
            self.run_regression()

            # c. error
            error = (sol.beta0-par.beta0_target)**2 + (sol.beta1-par.beta1_target)**2

            print(f'{par.alpha = :12.8f}, {par.sigma = :12.8f}: {error = :12.8f}')

            return error
        
        # b. initial guess
        if alpha is None: alpha = par.alpha
        if sigma is None: sigma = par.sigma

        x0 = np.array([alpha,sigma])
        results = optimize.minimize(obj,x0,method='Nelder-Mead',bounds=((0.001,0.999),(0.001,10)))

        assert results.success

        # c. final evaluation
        obj(results.x)