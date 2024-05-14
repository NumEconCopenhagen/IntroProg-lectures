
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
    
    def L_star_numerical(self,G=1,ufunc='cb'):
        par = self.par 
        sol = getattr(self,'sol_'+ufunc)

        # Define objective function
        u = getattr(self,'u_'+ufunc)

        obj = lambda L: -u(L,G)

        # Optimize
        sol.L_star = optimize.minimize_scalar(obj,bounds=(0,24)).x
          

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
    
    def utility_budget_cb(self):
        '''
        Compute utility at optimal labor supply given that the budget is balanced 
        '''
        par = self.par
        sol = self.sol_cb


        self.L_star_analytical()
        sol.G = par.tau*par.w*sol.L_analytical
        sol.U_star =self.u_cb(sol.L_analytical,sol.G)

    def tau_star_cb(self):
        '''
        Compute optimal tax rate for cobb douglas 
        '''
        par = self.par
        sol = self.sol_cb

        def obj(tau):
            par.tau = tau
            self.utility_budget_cb()
            return -sol.U_star
        
        sol.tau_star = optimize.minimize_scalar(obj,bounds=(0,1)).x
    

    def find_G_ces(self,tau=0.3,printit=True):
        '''
        Compute the G that balances the budget in the ces case
        '''
        par = self.par
        sol = self.sol_ces

        par.tau = tau 
        def obj(G):
            self.L_star_numerical(G,ufunc='ces')

            return G - par.tau * par.w * sol.L_star

        sol.G_star = optimize.root_scalar(obj,bracket=[0,tau*par.w*24  ]).root
        if printit:
            print(f'G_star = {sol.G_star:.2f}')

    def find_tau_star_ces(self):
        '''
        Compute optimal tax rate for CES
        '''
        par = self.par
        sol = self.sol_ces

        def obj(tau):
            par.tau = tau
            self.find_G_ces(tau,printit=False)
            U = self.u_ces(sol.L_star,sol.G_star)
            return -U
        
        sol.tau_star = optimize.minimize_scalar(obj,bounds=(0,1)).x
        print(f'tau_star = {sol.tau_star:.2f}')

    def plot_utiltiy(self,G=1,ufunc='cb'):
        '''
        Plot utility as a function of labor supply
        '''
        par = self.par
        sol = self.sol_cb

        # Create figure
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        # Create grid
        L_vec = np.linspace(10, 20, 100)

        # Compute utility
        u_vec = self.u_ces(L_vec, G)

        # Analytical solution
        self.L_star_analytical(G)

        # Numerical solution
        self.L_star_numerical(G,ufunc=ufunc)

        # Plot
        ax.plot(L_vec, u_vec)
        ax.axvline(sol.L_analytical, color='red', label='Analytical solution')
        ax.axvline(sol.L_star, color='green', linestyle='--', label='Numerical solution')
        ax.set_xlabel('Labor supply')
        ax.set_ylabel('Utility')
        ax.legend()
        ax.grid(True)
        plt.show()



    def plot_L_star_across_w(self,G=1):
        '''
        Plot optimal labor supply as a function of w
        '''
        par = self.par
        sol = self.sol_cb

        # Create figure
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        # Create grid
        w_vec = np.linspace(0.5, 3, 100)


        w_org = par.w
        # Compute optimal labor supply
        L_star_vec = np.empty_like(w_vec)
        for i, w in enumerate(w_vec):
            par.w = w
            self.L_star_analytical(G)
            L_star_vec[i] = self.sol_cb.L_analytical

        # Reset w
        par.w = w_org

        # Plot
        ax.plot(w_vec, L_star_vec)
        ax.set_xlabel('w')
        ax.set_ylabel('Labor supply')
        ax.grid(True)
        plt.show()


    def plot_model_across_tau_cb(self):
        '''
        Plot solution across 
        '''

        par = self.par
        sol = self.sol_cb

        tau_grid = np.linspace(0.1,0.9, 100)

        L_star_vec = np.empty_like(tau_grid)
        G_star_vec = np.empty_like(tau_grid)
        U_star_vec = np.empty_like(tau_grid)


        org_tau = par.tau
        for i, tau in enumerate(tau_grid):
            par.tau = tau
            self.utility_budget_cb()
            
            L_star_vec[i] = sol.L_analytical
            G_star_vec[i] = sol.G
            U_star_vec[i] = sol.U_star

        # reset tau
        par.tau = org_tau

        fig,axes = plt.subplots(3,1,figsize=(10,10))
        axes[0].plot(tau_grid,L_star_vec)
        axes[0].set_ylabel('Labor supply')
        axes[0].set_title('Labor supply')
        

        axes[1].plot(tau_grid,G_star_vec)
        axes[1].set_ylabel('Government spending')
        axes[1].set_title('Government spending')
        
        axes[2].plot(tau_grid,U_star_vec)
        axes[2].set_ylabel('Utility')
        axes[2].set_title('Utility')
        fig.tight_layout()
        plt.show()
        