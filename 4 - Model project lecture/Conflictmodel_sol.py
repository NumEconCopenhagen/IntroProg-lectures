
# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize, interpolate
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})




# Define the model
class ConflictModel():

    def __init__(self, **kwargs):
        '''
        Initialize the model with default parameters
        kwargs allow any parameter in the par namespace to be overridden by the user
        '''

        self.par = par = SimpleNamespace() # Create a namespace object for parameters
        self.sol = sol = SimpleNamespace() # Create a namespace object for solution results
        

        # Set default parameters
        self.setup()
        
        # Set default utility functions 
        self.u = self.u_ql

        # Update parameters with user input
        for key, value in kwargs.items():
            setattr(par, key, value)

        

    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par
        sol = self.sol

        # Model parameters
        par.ea = 10
        par.eb = 10
        par.epsilon = 10.
        par.eta = 1000. 

        par.small = 1e-4



    def u_ql(self,c,cm):
        '''
        Utility function
        '''
        par = self.par
        return c + cm**(1-1/par.epsilon)/(1-1/par.epsilon)
    
    def u_nonlinear(self,c,cm):
        par = self.par
        return c**(1-1/par.eta)/(1-1/par.eta) + cm**(1-1/par.epsilon)/(1-1/par.epsilon)
            


    def solve_A(self,p):
        par = self.par

        def obj(cm):
            return -self.u(par.ea-p*cm,cm)
        
        res = optimize.minimize_scalar(obj,bounds=(par.small ,par.ea/p-par.small),method='bounded')
        return res.x
    

    def solve_A_analytical(self,p):
        par = self.par
        return p**(-par.epsilon)



    def solve_A_grid(self):
        par = self.par 
        sol = self.sol

        sol.p_vec = np.linspace(1,2,500)
        sol.D_vec = np.empty_like(sol.p_vec)

        for i, p in enumerate(sol.p_vec):
            sol.D_vec[i] = self.solve_A(p)
        
        sol.D_func = interpolate.RegularGridInterpolator([sol.p_vec],sol.D_vec,method='linear')

        


    def solve_B(self,method='optimize'):
        par = self.par
        sol = self.sol

        if method in ['analytical']:
            D_func = lambda p : self.solve_A_analytical(p)  


        elif method in ['optimize']:
            D_func = lambda p : self.solve_A(p)  



        elif method in ['interpolate']:
            
            self.solve_A_grid()
            D_func = lambda p : sol.D_func([p]).item()

        def obj(p):
            Dp = D_func(p)
            c = np.fmax(par.eb-Dp,par.small)
            cm = np.fmax(p*Dp,par.small)

            return -self.u(c,cm)
        
        res = optimize.minimize_scalar(obj,bounds=(1,2),method='bounded')
        return res.x

    def solve_B_analytical(self):
        par = self.par
        return (par.epsilon/(par.epsilon-1))**(par.epsilon/(2*par.epsilon-1))



    def plot_solve_A(self):
        par = self.par
        
        p = np.linspace(1,2,100)

        cm = np.empty_like(p)
        cm_a = np.empty_like(p)

        for i,pi in enumerate(p):
            cm[i] = self.solve_A(pi)
            cm_a[i] = self.solve_A_analytical(pi)
        

        fig, ax = plt.subplots()
        ax.plot(p,cm,label='Numerical')
        ax.plot(p,cm_a,label='Analytical',linestyle='--')
        ax.legend()
        ax.set_xlabel('Price')
        ax.set_ylabel("Consumption of $c'$")
        ax.set_title("Optimal consumption of $c'$ for A as a function of price")
        plt.show()



    def plot_solve_B(self,xrange= [5,15]):
        par = self.par

        epsilon_org = par.epsilon


        epsilon_vec = np.linspace(*xrange,100)

        p_vec = np.empty_like(epsilon_vec)
        p_a_vec = np.empty_like(epsilon_vec)
        

        for i,epsilon in enumerate(epsilon_vec):
            par.epsilon = epsilon
            p_vec[i] = self.solve_B()
            p_a_vec[i] = self.solve_B_analytical()
        
        par.epsilon = epsilon_org

        fig, ax = plt.subplots()
        ax.plot(epsilon_vec,p_vec,label='Numerical')
        ax.plot(epsilon_vec,p_a_vec,label='Analytical',linestyle='--')
        ax.legend()
        ax.set_xlabel('Epsilon')
        ax.set_ylabel("Price")
        ax.set_title("Price as a function of epsilon for B")


    def plot_solve_B_range_ea(self):
        par = self.par

        ea_org = par.ea

        ea_vec = np.linspace(2.5,30,100)

        p_vec = np.empty_like(ea_vec)

        for i,ea in enumerate(ea_vec):
            par.ea = ea
            p_vec[i] = self.solve_B()
            

        par.ea = ea_org

        fig, ax = plt.subplots()
        ax.plot(ea_vec,p_vec,label='Numerical')
        ax.legend()
        ax.set_xlabel('Ea')
        ax.set_ylabel("Price")
        ax.set_title("Price as a function of ea for B in the range 2.5 to 30")
        plt.show()