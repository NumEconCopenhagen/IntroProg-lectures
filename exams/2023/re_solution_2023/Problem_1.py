from types import SimpleNamespace
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

class Consumer():

    def __init__(self,par):
        """ initialize """
    
        self.par = par
        self.sol = SimpleNamespace()
        self.sim = SimpleNamespace()

        self.allocate()

    def allocate(self):
        """ allocate memory """

        sim = self.sim  
        par = self.par

        sim.psi2 = np.empty(par.N)
        sim.nu2 = np.empty(par.N)
        sim.y2 = np.empty(par.N)

    def utility(self,c):
        """ utility function"""

        par = self.par
        return c**(1-par.rho) / (1-self.par.rho)

    def utility_period1(self,c1):
        """ utility in period 1 """

        par = self.par
        sim = self.sim

        a1 = par.m1-c1
        m2 = (1+par.r)*a1 + sim.y2

        return (self.utility(c1) + par.beta*np.mean(self.utility(m2)) ) 


    def simulate(self,seed= None):
        """ simulate the model """

        par = self.par
        sim = self.sim

        np.random.seed(seed)

        sim.psi2[:] = np.random.normal(0,1, par.N)
        sim.nu2[:] = np.random.uniform(low=0,high=1,size=par.N)

        self.calc_y2()

    def calc_y2(self):
        """ calculate y2 """

        par = self.par
        sim = self.sim

        sim.y2[:] = np.exp(-0.5*par.sigma_psi**2+par.sigma_psi*sim.psi2)
        I = sim.nu2<par.pi
        sim.y2[I] = 1.0


    def solve(self,print_sol=True):
        """ solve the model """

        par = self.par
        sim = self.sim
        sol = self.sol

        def obj(c):
            return -self.utility_period1(c)

        res = optimize.minimize_scalar(obj, method='bounded', bounds=(1e-8,par.m1))

        sol.c1 = res.x
        sol.v1 = -res.fun
        
        if print_sol:
            print(f'c1 = {sol.c1:6.3f}')
            print(f'v1 = {sol.v1:6.3f}')

    def constant_pars(self,plotit=True,seed=2020):
        """ find the combinations of pi and sigma_psi that gives the same expected utility in the first period """

        par = self.par
        sol = self.sol

        self.simulate(seed=seed)
        self.solve(print_sol=False)

        sol.pis = np.linspace(0.0,0.9,50)
        sol.sigma_psis = np.zeros(sol.pis.shape)

        v1_target = sol.v1
        pi_org = par.pi
        sigma_psi_org = par.sigma_psi        
        for i,pi in enumerate(sol.pis):

            def obj_equal(sigma_psi):
                
                par.sigma_psi = sigma_psi
                self.calc_y2()
                self.solve(print_sol=False)

                return sol.v1-v1_target

            par.pi = pi
            
            res = optimize.root_scalar(obj_equal, method='brentq', bracket=[0.0,1.00])
            sol.sigma_psis[i] = res.root  

        par.pi = pi_org
        par.sigma_psi = sigma_psi_org

        if plotit: self.plot_constant_pars()

    ################
    ### plotting ###
    ################

    def plot_y2(self):

        sim = self.sim
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(sim.y2, bins=50, density=True);

    def plot_v1(self):

        par = self.par
        sol = self.sol

        self.solve(print_sol=False)
        
        c1_vec = np.linspace(0.5,par.m1,100)
        v1_vec = [ self.utility_period1(c1) for c1 in c1_vec ]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(c1_vec,v1_vec,label='v1')
        ax.plot(sol.c1,sol.v1,'o',label='optimal choice')
        ax.set_xlabel('$c_1$')
        ax.set_ylabel('$v_1$')
        ax.legend()
        
    def plot_constant_pars(self):

        sol = self.sol

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(sol.pis,sol.sigma_psis)
        ax.set_xlabel('$\pi$')
        ax.set_ylabel('$\sigma_{\psi}$')
        ax.set_title('constant expected utility')
        ax.grid(True)

    def plot_v1_3dplot(self):

        par = self.par
        sol = self.sol

        # calculate values
        pi_org = par.pi
        sigma_psi_org = par.sigma_psi
       
        sol.v1_mat = np.empty((len(sol.pis),len(sol.sigma_psis)))
        for i,pi in enumerate(sol.pis):
            for j,sigma_psi in enumerate(sol.sigma_psis):

                par.pi = pi
                par.sigma_psi = sigma_psi

                self.calc_y2()
                self.solve(print_sol=False)

                sol.v1_mat[i,j] = sol.v1
                        
        par.pi = pi_org
        par.sigma_psi = sigma_psi_org

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        P, S = np.meshgrid(sol.pis,sol.sigma_psis)
        
        surf = ax.plot_surface(P, S, sol.v1_mat,label='Expected utility',alpha=0.7)
        surf._edgecolors2d = surf._edgecolor3d # For color in legend
        surf._facecolors2d = surf._facecolor3d # For color in legend
        
        ax.plot(sol.pis,sol.sigma_psis,np.diagonal(sol.v1_mat),color='red',linewidth=3,label='Constant expected utility')


        ax.set_xlabel('$\pi$')
        ax.set_ylabel('$\sigma_{\psi}$')
        ax.set_zlabel('$v_1$')
        ax.legend()
        ax.set_title('expected utility in period 1')

        fig.tight_layout()

        plt.show()