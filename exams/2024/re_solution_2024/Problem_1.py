
import numpy as np
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scipy import optimize


class DiscreteJobSearchModel:
    def __init__(self,**kwargs):
        '''
        Initialize the model with parameters
        '''
        self.par = par = SimpleNamespace()
        self.sol = sol = SimpleNamespace()


        par.max_iter = 1000

        for key,value in kwargs.items():
            setattr(par,key,value)
    

    def Ve(self,wk):
        par = self.par
        
        return self.u(wk) /(1-par.beta )

    def u(self,c):
        par = self.par
        return c**(1-par.rho)/(1-par.rho)

    
    def Vu_x(self,x,Vu_old=0):
        '''
        Vu for a given x
        '''
        par = self.par

        exp = par.pi[x:]@self.Ve(par.w[x:]) + np.sum(par.pi[:x])*Vu_old

        return self.u(par.z)+ par.beta * exp

    def Vu_max(self,Vu_old=0):
        '''
        Vu given that x is chosen optimally
        '''
        par = self.par

        exp = np.sum(par.pi*np.fmax(self.Ve(par.w),Vu_old))

        return self.u(par.z) + par.beta * exp

    def pi_func(self,e):
        par = self.par
        return np.exp(-par.w/(0.5+e))/np.sum(np.exp(-par.w/(0.5+e)))

       
    def Vu2_x(self,x,Vu_old=0,e=0):
        '''
        Vu in the search effort model for a given x and e
        '''
        par = self.par

        pi = self.pi_func(e)
    
        exp = pi[x:]@self.Ve(par.w[x:]) + np.sum(pi[:x])*Vu_old

        return self.u(par.z - e) + par.beta * exp
    

    def Vu2_max(self,Vu_old=0,e=0):
        '''
        Vu in the search effort model, given that x is chosen optimally
        '''
        par = self.par

        pi = self.pi_func(e)

        exp = np.sum(pi*np.fmax(self.Ve(par.w),Vu_old))

        return  self.u(par.z - e)+ par.beta * exp
    
    def find_x(self,Vu_next,printit=False):
        '''
        Find the optimal x given Vu_next
        '''
        par = self.par
        sol = self.sol


        Vu_vec = np.array([self.Vu_x(x,Vu_next) for x in range(par.K+1)])
        x_new = np.argmax(Vu_vec)
        Vu = Vu_vec[x_new]

        if printit:
            
            print(f'Vu = {Vu:.4f}')
            print(f'x =  {x_new+1}') # Python indexing starts a 0.
        else:
            return x_new, Vu


    def find_x_e(self,Vu_new,printit=False):
        '''
        Find the optimal x and e given Vu in the search effort model
        '''

        par = self.par
        sol = self.sol


        def obj(e):
            Vu_vec = np.array([self.Vu2_x(x,Vu_new,e) for x in range(par.K+1)])
        
            return -np.max(Vu_vec)

        opt_sol = optimize.minimize_scalar(obj,bounds=(0,par.z),method='bounded') 
        
        e_new = opt_sol.x
        Vu = -opt_sol.fun
        x_new = np.argmax(np.array([self.Vu2_x(x,Vu_new,e_new) for x in range(par.K+1)]))
        if printit:
            
            print(f'Vu = {Vu:.4f}')
            print(f'x =  {x_new+1}')
            print(f'e =  {e_new:.4f}')
        else:
            return x_new,e_new, Vu


    def find_e(self,Vu_new):
        '''
        Find the optimal value of e, given that x is chosen optimally in the search effort model
        '''
        par = self.par
        sol = self.sol


        def obj(e):
        
            return -self.Vu2_max(Vu_new,e)

        opt_sol = optimize.minimize_scalar(obj,bounds=(0,par.z),method='bounded') 
        
        e_new = opt_sol.x
        Vu = -opt_sol.fun

        return e_new, Vu
                

    def solve_bellman(self, print_iter=False,model2=False,Vu_guess=0):
        '''
        Solve the bellman equation to find Vu*
        '''

        par = self.par
        sol = self.sol


        Vu_new = Vu_guess
        # Solve
        for i in range(par.max_iter):
            
            if model2:
                e_new, Vu = self.find_e(Vu_new)

            else:
                Vu =self.Vu_max(Vu_new)
                
            
            if print_iter:
                if model2:
                    print(f'Iter {i:3d}:  Vu = {Vu:.4f}, e = {e_new:.4f}')    
                else:
                    print(f'Iter {i:3d}:  Vu = {Vu:.4f}')

            if np.abs(Vu_new-Vu)<par.tol:
                break
            
            Vu_new = Vu

        if i == par.max_iter-1:
            print('No convergence')
            sol.solved = False
            return
        else:

            sol.V = Vu_new
            sol.x = self.find_x_e(Vu_new)[0]
            sol.solved = True
            sol.iter = i

            if model2:
                sol.e = e_new
    

    def plot_sol_across_z(self,N=100,model2=False):
        '''
        Plot solution values of 
        '''

        par = self.par
        sol = self.sol


        z_vec = np.linspace(0,4,N)

        xstar_vec = np.empty(N)
        Vustar_vec = np.empty(N)
        S_vec = np.empty(N)

        z_org = par.z
        for i,z in enumerate(z_vec):
            par.z = z
            self.solve_bellman()
            xstar_vec[i] = sol.x +1 
            Vustar_vec[i] = sol.V
            S_vec[i] = par.beta * (par.pi[sol.x:]@self.Ve(par.w[sol.x:]) )

        if model2:
            fig,ax = plt.subplots(2,2,figsize=(10,10),sharex=True)
            ax1 = ax[0,0]
            ax2 = ax[1,0]
            ax3 = ax[0,1]
            ax4 = ax[1,1]

            ax4.plot(z_vec,[0]*N)
            ax4.set_title(r'$e^*$')
            ax4.set_xlabel('z')
            ax2.set_xlabel('z')
        
        else:
            fig,ax = plt.subplots(3,1,figsize=(6,10),sharex=True)
            ax1 = ax[0]
            ax2 = ax[1]
            ax3 = ax[2]
            ax3.set_xlabel('z')

        ax1.plot(z_vec,xstar_vec)
        ax1.set_title(r'$x^*$')
        ax2.plot(z_vec,Vustar_vec)
        ax2.set_title(r'$V_u^*$')
        ax3.plot(z_vec,S_vec)
        ax3.set_title('S')
        
            
        if model2:
            e_vec = np.empty(N)
            for i,z in enumerate(z_vec):
                par.z = z
                self.solve_bellman(model2=True)
                xstar_vec[i] = sol.x +1
                Vustar_vec[i] = sol.V
                e_vec[i] = sol.e
                pi = self.pi_func(sol.e)
                S_vec[i] = par.beta * (pi[sol.x:]@self.Ve(par.w[sol.x:]) )
                
            
            ax1.plot(z_vec,xstar_vec)
            ax2.plot(z_vec,Vustar_vec)
            ax3.plot(z_vec,S_vec)
            ax4.plot(z_vec,e_vec)

            fig.legend(['Original model','Search effort model'])
        
        
        fig.tight_layout()
        par.z = z_org


    def plot_sol_across_fixede(self,N=100):
        par = self.par
        sol = self.sol


        e_vec = np.linspace(0.1,2,N)

        xstar_vec = np.empty(N)
        Vustar_vec = np.empty(N)
        S_vec = np.empty(N)

        pi_org = par.pi
        for i,e in enumerate(e_vec):
            par.pi = np.exp(-par.w/e)/np.sum(np.exp(-par.w/e))
            self.solve_bellman()
            xstar_vec[i] = sol.x +1
            Vustar_vec[i] = sol.V
            S_vec[i] = par.beta * (par.pi[sol.x:]@self.Ve(par.w[sol.x:]) )

        par.pi = pi_org

        fig,ax = plt.subplots(3,1,figsize=(10,8))
        ax[0].plot(e_vec,xstar_vec)
        ax[0].set_title(r'$x^*$')
        ax[1].plot(e_vec,Vustar_vec)
        ax[1].set_title(r'$V_u^*$')
        ax[2].plot(e_vec,S_vec)
        ax[2].set_title('S')

        fig.tight_layout()




#### Solution functions that explictely solve for x in each iteration
        
    def solve_bellman_old(self, print_iter=False,model2=False):
        par = self.par
        sol = self.sol



        Vu_new = 1
        # Solve
        for i in range(par.max_iter):
            

            if model2:
                x_new, e_new, Vu = self.find_x_e(Vu_new)
            else:
                x_new, Vu = self.find_x(Vu_new)
            
            
            if print_iter:
                if model2:
                    print(f'Iter {i:3d}:  Vu = {Vu:.4f}, x = {x_new+1:3d}, e = {e_new:.4f}')    
                else:
                    print(f'Iter {i:3d}:  Vu = {Vu:.4f}, x = {x_new+1}')

            if np.isclose(Vu_new,Vu):
                break
            
            Vu_new = Vu

        if i == par.max_iter-1:
            print('No convergence')
            sol.solved = False
            return
        else:
            sol.V = Vu_new
            sol.x = x_new
            sol.solved = True
            sol.iter = i

            if model2:
                sol.e = e_new