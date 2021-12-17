from types import SimpleNamespace
import numpy as np
from scipy import optimize

class WorkerCapitalistEconomyClass():

    def __init__(self):

        par = self.par = SimpleNamespace()

        # a. parameters
        par.kappa = 0.1 # home production
        par.omega = 10 # disutility of labor supply factor
        par.eta = 1.50 # curvature of disutility of labor supply
        par.alpha = 0.50 # curvature of production function
        par.Nw = 99 # number of workers
        par.Nc = 1 # number of capitalists

        # b. grids
        par.num_w = 10
        par.grid_w = np.linspace(0.1,1.5,par.num_w)
        par.grid_mkt_clearing = np.zeros(par.num_w)

        # c. initial prices
        self.p = 1 # output price
        self.w = 1 # wage

    def utility_w(self,c,l):
        """ utility of workers """
        
        par = self.par

        return np.log(c+par.kappa)-par.omega*l**par.eta

    def workers(self):
        """ maximize utility for workers """
        
        p = self.p
        w = self.w

        # a. solve
        obj = lambda l: -self.utility_w((w*l)/p,l) # substitute in the budget constraint
        res = optimize.minimize_scalar(obj,bounds=(0,1),method='bounded')
        
        # b. save
        self.l_w_star = res.x
        self.c_w_star = (w*self.l_w_star)/p
        
    def utility_c(self,c,l):
        """ utility of capitalists """

        par = self.par

        return np.log(c+par.kappa)-par.omega*l**par.eta

    def capitalists(self):
        """ maximize utility of capitalists """

        p = self.p
        w = self.w
        pi = self.pi

        # a. solve
        obj = lambda l: -self.utility_c((w*l+pi)/p,l) # subsittute in the budget constraint
        res = optimize.minimize_scalar(obj,bounds=(0,1),method='bounded')
        
        # b. save
        self.l_c_star = res.x
        self.c_c_star = (w*self.l_c_star+pi)/p
        
    def firm(self):
        """ maximize firm profits """
        
        p = self.p
        w = self.w
        par = self.par

        # a. solve
        f = lambda l: l**par.alpha
        obj = lambda l: -(p*f(l)-w*l)
        x0 = [0.0]
        res = optimize.minimize(obj,x0,bounds=((0,None),),method='L-BFGS-B')
        
        # b. save
        self.l_star = res.x[0]
        self.y_star = f(self.l_star)
        self.Pi = p*self.y_star - w*self.l_star
        
    def evaluate_equilibrium(self):
        """ evaluate equilirium """

        par = self.par

        # a. optimal behavior of firm
        self.firm()
        self.pi = self.Pi/par.Nc

        # b. optimal behavior of households
        self.workers()
        self.capitalists()

        # c. market clearing
        self.goods_mkt_clearing = par.Nw*self.c_w_star + par.Nc*self.c_c_star - self.y_star
        self.labor_mkt_clearing = par.Nw*self.l_w_star + par.Nc*self.l_c_star - self.l_star
    
    def find_equilibrium(self):

        par = self.par

        # a. grid search
        print('grid search:')
        for i,w in enumerate(par.grid_w):
            self.w = w
            self.evaluate_equilibrium()
            par.grid_mkt_clearing[i] = self.goods_mkt_clearing
            print(f' w = {w:.2f} -> {par.grid_mkt_clearing[i]:12.8f}')
        
        print('')

        # b. find bounds
        left = np.max(par.grid_w[par.grid_mkt_clearing < 0])
        right = np.min(par.grid_w[par.grid_mkt_clearing > 0])
        print(f'equilibrium price must be in [{left:.2f},{right:.2f}]\n')            

        # c. bisection search
        def obj(w):
            self.w = w
            self.evaluate_equilibrium()
            return self.goods_mkt_clearing

        res = optimize.root_scalar(obj,bracket=[left,right],method='bisect')
        self.w_eq = res.root
        print(f'the equilibrium wage is {self.w_eq:.4f}\n')

        # d. show result
        u_w = self.utility_w(self.c_w_star,self.l_w_star)
        print(f'workers      : c = {self.c_w_star:6.4f}, l = {self.l_w_star:6.4f}, u = {u_w:7.4f}')
        u_c = self.utility_c(self.c_c_star,self.l_c_star)
        print(f'capitalists  : c = {self.c_c_star:6.4f}, l = {self.l_c_star:6.4f}, u = {u_c:7.4f}')        
        print(f'goods market : {self.goods_mkt_clearing:.8f}')
        print(f'labor market : {self.labor_mkt_clearing:.8f}')