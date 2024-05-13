
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



    def u_cb(self, L,G):
        '''
        Utility function
        '''
        par = self.par
        C = par.kappa + (1-par.tau)*par.w*L
        C_bounded = np.fmax(C,par.zero_tol)
        G_bounded = np.fmax(G,par.zero_tol)

        return np.log(C_bounded**par.alpha * G_bounded**(1-par.alpha)) - par.nu*L**2/2


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
    

    def solve_L(self,G=1,returnit=False,utype='cb'):
        '''
        Solve for optimal labor supply
        '''
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')


        u_func = getattr(self,f'u_{utype}')
        # Define objective function
        obj = lambda L: -u_func(L,G)

        # Call optimizer
        sol.L = optimize.minimize_scalar(obj,bounds=[0,24]).x

        if returnit:
            return sol.L

    def solve_L_analytical(self,G=1,returnit=False):
        '''
        Solve for optimal labor supply analytically
        '''
        par = self.par
        sol_cb = self.sol_cb
        

        wtilde = (1-par.tau)*par.w

        sol_cb.L_analytical =(-par.kappa+ np.sqrt(par.kappa**2 + (4*par.alpha/par.nu)*wtilde**2))/(2*wtilde)

        if returnit:
            return sol_cb.L_analytical
        

    def Q1(self):

        par = self.par
        sol_cb = self.sol_cb

        for G in [1.,2.]:
            # Solve for optimal labor supply
            self.solve_L(G=G)
            # Solve for optimal labor supply analytically
            self.solve_L_analytical(G=G)

            # Print results
            print(f'G = {G}:')
            print(f'Numerical solution: L = {sol_cb.L:.4f}')
            print(f'Analytical solution: L = {sol_cb.L_analytical:.4f}')



    def Q3(self,tau_bounds = [0,1],n=100):

        par = self.par
        sol_cb = self.sol_cb

        # Define a grid of tau values
        tau_values = np.linspace(*tau_bounds,n)

        # Initialize arrays to store L, G, and utility values
        L_values = np.empty_like(tau_values)
        G_values = np.empty_like(tau_values)
        utility_values = np.empty_like(tau_values)

        tau_org =par.tau
        # For each tau value, calculate L, G, and utility
        for i, tau in enumerate(tau_values):

            par.tau = tau
            L_values[i] = self.solve_L(G=1, returnit=True)
            G_values[i] = par.tau * par.w * L_values[i]
            utility_values[i] = self.u_cb(L_values[i], G_values[i])
        
        par.tau = tau_org

        # Create subplots
        fig, axes = plt.subplots(3, 1, figsize=(8, 12))

        # Plot L, G, and utility against tau
        axes[0].plot(tau_values, L_values, label='L')
        axes[0].set_title('L vs tau')
        axes[0].set_xlabel('tau')
        axes[0].set_ylabel('L')

        axes[1].plot(tau_values, G_values, label='G')
        axes[1].set_title('G vs tau')
        axes[1].set_xlabel('tau')
        axes[1].set_ylabel('G')

        axes[2].plot(tau_values, utility_values, label='Utility')
        axes[2].set_title('Utility vs tau')
        axes[2].set_xlabel('tau')
        axes[2].set_ylabel('Utility')

        plt.tight_layout()
        plt.show()

        

    def find_tau(self,G=1,utype='cb'):
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')
        u_func = getattr(self,f'u_{utype}')

        if utype in ['cb']:
            def obj(tau):
                par.tau = tau
                self.solve_L(G=G,utype=utype)
                u = u_func(sol.L,par.tau * par.w * sol.L)
                return -u # Minimize negative utility
    
        elif utype in ['ces']:
            def obj(tau):
                par.tau =tau 
                G = self.find_G(plotit=False,returnit=True,printit=False,utype=utype)
                self.solve_L(G=G,utype=utype)
                u = u_func(sol.L, par.tau * par.w * sol.L)
                return -u # Minimize negative utility
        par.tau = optimize.minimize_scalar(obj,bounds=[0,1]).x
        

        # resolve L
        G = self.find_G(plotit=False,returnit=True,printit=False,utype=utype)
        self.solve_L(G=G,utype=utype)
        print(f'The optimal tax rate is {par.tau:.4f}')
        print(f'The optimal labor supply is {sol.L:.4f}')
        G_out = par.tau * par.w * sol.L
        assert np.isclose(G,G_out), 'The budget is not balanced'
        print(f'The optimal government spending is {G:.4f}')
        print(f'The optimal utility is {u_func(sol.L,G):.4f}')

        sol.G_opt = G
        sol.tau_opt = par.tau
        return par.tau


    def find_G(self,plotit=True,printit=True,returnit=False,utype='ces'):
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')

        if utype in ['cb']:
            sol.G_opt = par.tau * par.w * sol.L
            if returnit:
                return sol.G_opt  
            else:
                return None 
            
        def obj(G):
            
            self.solve_L(G=G,utype=utype)
            G_next = par.tau * par.w * sol.L

            return G-G_next

        G_opt = optimize.root_scalar(obj,x0= sol_cb.G_opt, x1= sol_cb.G_opt*0.9).root
        if printit:
            print(f'{G_opt:.4f} is the budget neutral government spending')


        if plotit:
            # Plot
            G_grid = np.linspace(0,10,100)
            L_grid = np.empty_like(G_grid)
            budget_grid = np.empty_like(G_grid)

            for i,G in enumerate(G_grid):
                self.solve_L(G=G,utype=utype)
                L_grid[i] = sol.L
                budget_grid[i] = G - (par.tau)*par.w*sol.L
            
            fig,ax = plt.subplots(2,1,figsize=(6,12))
            ax[0].plot(G_grid,L_grid)
            ax[0].set_xlabel('G')
            ax[0].set_ylabel('L')
            ax[0].set_title('Labor supply as a function of G')

            ax[1].plot(G_grid,budget_grid)
            ax[1].set_xlabel('G')
            ax[1].set_ylabel('Budget surplus')
            ax[1].set_title('Budget surplus as a function of G')

            ax[0].axvline(G_opt,color='red',linestyle='--')
            ax[1].axvline(G_opt,color='red',linestyle='--')

        if returnit:
            return G_opt    





    def verify_L_solve(self,G=1,L_bounds=[10,20],utype='cb'):
        '''
        Verify that the numerical solution is correct
        '''
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')
        L_vec = np.linspace(*L_bounds,100)
        u_func = getattr(self,f'u_{utype}')

        # Calculate utility for each L
        u_vec = np.empty_like(L_vec)
        for i,L in enumerate(L_vec):
            u_vec[i] = u_func(L,G)

        # Find maximum
        self.solve_L(G,utype=utype)

        # Plot
        fig,ax = plt.subplots()
        ax.plot(L_vec,u_vec)
        ax.axvline(sol.L,color='red',linestyle='--')
        ax.set_xlabel('L')
        ax.set_ylabel('u')
        ax.set_title('Utility as a function of L')

    def plot_L(self,utype='cb',G=1):
        '''
        Plot labor supply as a function of w
        '''
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')

        w_org = par.w
        # Create grid
        w_vec = np.linspace(0.5,2,100)
        L_vec = np.empty_like(w_vec)
        for w in w_vec:
            par.w = w
            self.solve_L(G,utype=utype)
            L_vec[w_vec==w] = sol.L
        
        par.w = w_org

        # Plot
        fig,ax = plt.subplots()
        ax.plot(w_vec,L_vec)
        ax.set_xlabel('w')
        ax.set_ylabel(r'$L^*$')
        ax.set_title('Labor supply as a function of w')


    def plot_U_across_tau(self,utype='cb'):
        '''
        Plot utility as a function of tau
        '''
        par = self.par
        sol_cb = self.sol_cb
        sol_ces = self.sol_ces
        sol = getattr(self,f'sol_{utype}')
        u_func = getattr(self,f'u_{utype}')
        
        
        
        tau_org = par.tau
        # Create grid
        tau_vec = np.linspace(0.1,0.9,100)
        u_vec = np.empty_like(tau_vec)

        for i,tau in enumerate(tau_vec):
            par.tau = tau
            self.find_G(plotit=False,printit=False,utype=utype)
            u_vec[i] = u_func(sol.L,par.tau * par.w * sol.L)
        
        par.tau = tau_org

        # Plot
        fig,ax = plt.subplots()
        ax.plot(tau_vec,u_vec)
        ax.set_xlabel(r'$\tau$')
        ax.set_ylabel('u')
        ax.set_title('Utility as a function of tau')