# Import packages (just the ones you need)
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy import optimize
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})
from scipy import interpolate
from matplotlib import cm


# Define the model
class ConsumerModelClass2():

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
        self.draw_shocks()


    def setup(self):
        '''
        Set default parameters
        '''
        par = self.par

        # Model parameters
        par.alpha = 0.1
        par.rho = 2
        par.kappa = 0.5
        par.nu = 0.1
        par.r = 0.04
        par.beta = 0.94
        par.Delta = 0.5

        
        # Helpful variables 
        par.approx0 = 1e-5

        par.m2N = 200 
        par.m1N = 100
        par.d2N = 100

        par.simN = 10000


    def allocate(self):
        '''
        Allocate arrays for simulation and solution
        '''

        par = self.par
        sol = self.sol
        sim = self.sim

        # solution vecs
        sol.m1_vec = np.linspace(par.approx0*100,4,par.m1N)
        for vec in ['c1','v1','d1']:
            setattr(sol,f'{vec}_vec',np.empty(par.m1N))
        

        sol.m2_vec = np.linspace(1-par.Delta,5,par.m2N) # Set minimum bound for m2 1-Delta because that is lowest m2 with no savings 
        sol.d2_vec = np.linspace(par.approx0,5,par.d2N)


        for vec in ['c2','v2']:
            setattr(sol,f'{vec}_grid',np.empty((par.m2N,par.d2N)))


        # Simulation varibles 
        simvars = ['m1','m2','y2','c1','c2','d1','v1','v2']
        for var in simvars:
            setattr(sim,var,np.empty(par.simN))
        


    def draw_shocks(self,seed=1997):
        par = self.par
        sim = self.sim

        np.random.seed(seed)
        sim.m1[:] = np.fmax(np.random.normal(2,0.4,size=par.simN),0.1) #Set minimum bound for m1 0.1
        sim.y2[:] = np.random.choice([1-par.Delta,1+par.Delta],p=[0.5,0.5],size=par.simN)
            
    def utility(self, c,d):
        par = self.par
        return c**(1-par.rho)/(1-par.rho)+ par.alpha*d**(1-par.rho)/(1-par.rho)

    def bequest(self, m, c,d):
        par = self.par
        return par.nu*(m+d-c+par.kappa)**(1-par.rho)/(1-par.rho)

    def v2(self, c2, d2, m2):
        return self.utility(c2,d2) + self.bequest(m2, c2,d2)

    def v1(self, c1, d1, m1):
        par = self.par
        sol = self.sol 

        # a. v2 value, if low income
        m2_low = (1+par.r)*(m1-c1-d1) + 1-par.Delta
        v2_low = sol.v2_interp([m2_low,d1])[0]
        
        # b. v2 value, if high income
        m2_high = (1+par.r)*(m1-c1-d1) + 1+par.Delta
        v2_high = sol.v2_interp([m2_high,d1])[0]
        
        # c. expected v2 value
        v2 = 0.5*v2_low + 0.5*v2_high
        
        # d. total value
        return self.utility(c1,d1) + par.beta*v2
    
    def solve_period_2(self):
        par = self.par
        sol = self.sol        

        # b. solve for each m2 in grid
        for i,m2 in enumerate(sol.m2_vec):
            for j,d2 in enumerate(sol.d2_vec):

                # i. objective
                obj = lambda x: -self.v2(x[0],d2,m2)

                # ii. initial value (consume half)
                x0 = m2/2

                # iii. optimizer
                result = optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((par.approx0*1e-8,m2-par.approx0*1e-8),))

                # iv. save
                sol.v2_grid[i,j] = -result.fun
                sol.c2_grid[i,j] = result.x[0]
                
        

    def solve_period_1(self):
        par = self.par
        sol = self.sol

        # b. solve for each m1 in grid
        for i,m1 in enumerate(sol.m1_vec):
            
            # i. objective
            obj = lambda x: -self.v1(x[0],x[1],m1)

            # ii. initial guess (consume half)
            x0 = [m1/3,m1/3]
            
            # iii. optimize
            bound = (par.approx0*1e-4,m1-par.approx0*1e-4)
            bounds = (bound, bound)
            ineq_con = {'type': 'ineq', 'fun': lambda x: m1-x[0]-x[1]} 


            # This step can give the warning that x steps out of bounds, this is not a issue for us
            result = optimize.minimize(obj,x0, method='SLSQP',
                                    bounds=bounds,
                                    constraints=[ineq_con])
            
            
            # iv. save
            sol.v1_vec[i] = -result.fun
            sol.c1_vec[i] = result.x[0]
            sol.d1_vec[i] = result.x[1]
    

    def solve(self):
        par = self.par
        sol = self.sol

        # a. solve period 2
        self.solve_period_2()

        # b. construct interpolator
        sol.v2_interp = interpolate.RegularGridInterpolator((sol.m2_vec,sol.d2_vec), sol.v2_grid,
                                                    bounds_error=False,fill_value=None)
        # c. solve period 1
        self.solve_period_1()


    def plot_c1(self):
        par = self.par
        sol = self.sol

        fig,ax  = plt.subplots()
        ax.plot(sol.m1_vec,sol.c1_vec,label='non-durable consumption')
        ax.plot(sol.m1_vec,sol.d1_vec,label='durable consumption')
        ax.set_xlabel('$m_1$')
        ax.set_xlim([0,4])
        ax.set_ylim([0,2.5])
        ax.set_title('Consumption in period 1')
        ax.legend()

        plt.show()

    def plot_c2(self):
        par = self.par
        sol = self.sol


        m2_grid,d2_grid = np.meshgrid(sol.m2_vec,sol.d2_vec,indexing='ij')

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1,projection='3d')
        cs = ax.plot_surface(m2_grid,d2_grid,sol.c2_grid,cmap=cm.jet)

        # d. add labels
        ax.set_xlabel('$m_2$')
        ax.set_ylabel('$d_2$')


        # e. invert xaxis
        ax.invert_xaxis()

        ax.set_title('$c_2$')

        # f. add colorbar
        fig.colorbar(cs);



    def simulate(self):
        par = self.par
        sim = self.sim
        sol = self.sol

        # Interpolate solution
        sol.c2_interp = interpolate.RegularGridInterpolator((sol.m2_vec,sol.d2_vec), sol.c2_grid,
                                                    bounds_error=False,fill_value=None)
        sol.c1_interp = interpolate.RegularGridInterpolator([sol.m1_vec],sol.c1_vec,
                                                            bounds_error=False,fill_value=None)
        sol.d1_interp = interpolate.RegularGridInterpolator([sol.m1_vec],sol.d1_vec,
                                                            bounds_error=False,fill_value=None)
        sol.v1_interp = interpolate.RegularGridInterpolator([sol.m1_vec],sol.v1_vec,
                                                            bounds_error=False,fill_value=None)
        

        # Period 1
        sim.c1[:] = sol.c1_interp(sim.m1)
        sim.d1[:] = sol.d1_interp(sim.m1)
        sim.v1[:] = sol.v1_interp(sim.m1)

        # Period 2 
        sim.m2[:] = (1+par.r)*(sim.m1-sim.c1-sim.d1) + sim.y2
        sim.c2[:] = sol.c2_interp((sim.m2,sim.d1))
        sim.v2[:] = self.v2(sim.c2,sim.d1,sim.m2)

    def plot_simulation(self,axes = None,label=None,alpha=0.7):
        par = self.par 
        sim = self.sim

        if axes is None:
            fig,axes  = plt.subplots(figsize=(12,10),nrows = 3, ncols = 2)
        else:
            fig = axes[0,0].get_figure()
        
        # Prepare for making the plots with loops by making lists
        
        vecs = [sim.c1,sim.d1,sim.c2,sim.m1,sim.m2,sim.v1]
        names = [r'$c_{1}$',r'$d_{1}$',r'$c_{2}$',r'$m_{1}$',r'$m_{2}$',r'$v_{1}$']
        bins_list = [30,30,50,50,50,50]

        #plot using lists
        for ax, vec, name,bins in zip(axes.flatten(),vecs,names,bins_list):
            ax.hist(vec,bins=bins,density=True,label=label,alpha=alpha)
            ax.set_xlabel(name)
            ax.set_ylabel('Density')

        fig.tight_layout()

        return axes
