from types import SimpleNamespace
from scipy import optimize, interpolate
import numpy as np
import matplotlib.pyplot as plt

class ConflictModel:

    def __init__(self,par_dict={}) -> None:

        self.par = par = SimpleNamespace()
        self.sol = sol = SimpleNamespace()
        
        # parameters
        par.epsilon = 10.
        par.eta = 1000.
        par.sigma = 1.
        par.eA = 10.
        par.eB = 10.

        # interpolation options
        par.interp_N = 750
        par.p_buff = 1e-2
        par.small_number = 1e-4
        par.interp_range = [1-par.p_buff,1.6]

        # solution 
        sol.solved = False

        # set utliity function
        self.utility = self.utility_ql

        for key,value in par_dict.items():
            setattr(par,key,value)

    def __str__(self):
        
        text = ''
        for key,value in self.par.__dict__.items():
            text+= f'{key} = {value}\n'

        sol = self.sol
        if sol.solved:
            text += '\n' + self.print_solution()
            
        return text

    def __repr__(self):

        return self.__str__()

    def print_solution(self):
        
        sol = self.sol
        text = 'solution: \n'

        for key,value in sol.__dict__.items():

            if type(value) in [float]:
                text+= f'{key:9} = {value:3.3f}\n'
            elif type(value) in [int]:
                text+= f'{key:9} = {value}\n'

        return text


    def utility_ql(self,c,cm):

        par = self.par
        return c + (cm**(1-1/par.epsilon))/(1-1/par.epsilon)
    
    def utility_nonl(self,c,cm):

        par = self.par
        return  c**(1-1/par.eta)/(1-1/par.eta) + cm**(1-1/par.epsilon)/(1-1/par.epsilon)
        
    def analytical_quasilinear(self):
        """analytical solution for quasilinear utility function"""

        par = self.par
        p = (par.epsilon/(par.epsilon-1))**(par.epsilon/(2*par.epsilon-1))
        cm = p**(-par.epsilon)
        c = par.eA - p*cm

        print(f'p  = {p:.3f}')
        print('A consumes:')
        print(f"cm = {cm:.3f}")
        print(f'c  = {c:.3f}')

        c = par.eB - cm
        cm = p*cm 
        print('B consumes:')
        print(f"cm = {cm:.3f}")
        print(f'c  = {c:.3f}')
        
    def solve_buyer(self,p):

        par = self.par
        
        def obj(cm):
            c = np.fmax(par.eA - p*cm,par.small_number)
            return -self.utility(c,cm)
        
        res = optimize.minimize_scalar(obj,bounds=(par.small_number,par.eA/p),method='bounded')

        cm = res.x.item()
        if np.isclose(cm,par.small_number) or np.isclose(cm,par.eA/p):
            print('bounded solution')
            print(cm)
            print(p)

        return cm

    def solve_buyer_grid(self):

        par = self.par
        sol = self.sol

        # create demand function using interpolation
        sol.p_vec = np.linspace(*par.interp_range,par.interp_N)
        sol.D_vec = np.empty(par.interp_N)

        for i, p in enumerate(sol.p_vec):
            sol.D_vec[i] = self.solve_buyer(p,)

        sol.D_func = interpolate.RegularGridInterpolator([sol.p_vec],sol.D_vec,
                                                  bounds_error=False,
                                                  fill_value=None)

    def solve_seller(self,method='nested'):
        """solve seller problem using nested optimization or interpolation"""
        
        par = self.par
        sol = self.sol
        
        def obj(p):

            if method == 'interp':
            
                D_p = sol.D_func([p])
            
            elif method == 'nested':

                # Calculate demand for given price
                D_p = self.solve_buyer(p)
            
            c = par.eB -D_p
            cm = p*D_p

            # restrict p to positive consumption of both goods

            if c < par.small_number or cm < par.small_number:

                # if this is not held set both close to zero
                c = par.small_number
                cm = par.small_number

            return -self.utility(c,cm)

        p = optimize.minimize_scalar(obj,bounds=par.interp_range,method='bounded').x
        sol.p = p.item()

        if np.isclose(p,par.interp_range[0]) or np.isclose(p,par.interp_range[1]):
            print('bounded solution')
            print(p)
    
    def solve(self,print_sol=True,method='nested'):
        """
        solve model using interpolation over buyers demand function
        or nested optimization
        """

        par = self.par
        sol = self.sol

        if method == 'interp':

            # precalculate buyers demand function
            self.solve_buyer_grid()
        
        self.solve_seller(method=method)

        # buyer consumption
        sol.cm_buyer = self.solve_buyer(sol.p,)
        sol.c_buyer  = par.eA - sol.p*sol.cm_buyer
        
        # seller consumption
        sol.cm_seller = sol.p*sol.cm_buyer
        sol.c_seller = par.eB - sol.cm_buyer

        assert np.isclose(sol.c_buyer+sol.cm_seller,par.eA) and np.isclose(sol.c_seller+sol.cm_buyer,par.eB),'All endownments where not consumed'

        sol.solved = True

        if print_sol:
            print(self.print_solution())

    ############
    # plotting #
    ############

    def plot_demand(self):

        par = self.par
        sol = self.sol

        self.solve_buyer_grid()

        fig,ax = plt.subplots()
        ax.plot(sol.p_vec,sol.D_vec,label='numerical solution')
        
        if self.utility==self.utility_ql:

            # show analytical solution
            D_vec = sol.p_vec**(-par.epsilon)

            ax.plot(sol.p_vec,D_vec,linestyle='--',label='analytical solution')
        
        ax.legend(loc='upper right')
        ax.set_xlabel('p')
        ax.set_ylabel('D(p)')
        ax.set_title('demand function')

        plt.show()

    def plot_seller_utility(self):

        par = self.par
        sol = self.sol
        
        self.solve(method='interp')

        seller_utility_vec = self.utility(par.eB-sol.D_vec,sol.p_vec*sol.D_vec)
        
        fig,ax = plt.subplots()
        ax.plot(sol.p_vec,seller_utility_vec,label='seller utility')
        ax.vlines(sol.p,seller_utility_vec.min(),seller_utility_vec.max(),linestyle='--',label="seller's price",color='red')

        ax.legend(loc='lower right')
        ax.set_xlabel('p')
        ax.set_ylabel('U(p)')
        ax.set_title('seller utility')

        plt.show()

    def plot_p_across_epsilon(self):

        par = self.par
        sol = self.sol

        epsilon_vec = np.linspace(5,30,50)
        p_vec = np.empty_like(epsilon_vec)

        org_eps = par.epsilon
        for i, epsilon in enumerate(epsilon_vec):
            self.par.epsilon = epsilon
            self.solve(print_sol=False)
            p_vec[i] = self.sol.p

        par.epsilon = org_eps # Reset epsilon

        fig,ax = plt.subplots()
        
        ax.plot(epsilon_vec, p_vec,label='p')
        
        ax.set_xlabel(r'$\epsilon$')
        ax.set_ylabel('p*')
        ax.set_title('price')

        plt.show()

    def plot_p_across_eA(self):

        par = self.par
        sol = self.sol

        e_vec = np.linspace(2.5,30,25)
        p_vec_b = np.empty_like(e_vec)

        e_org = par.eA
        for i, e in enumerate(e_vec):

            # Sslve for e varying for buyer
            par.eA = e
            self.solve(print_sol=False)
            p_vec_b[i] = sol.p

        # reset 
        par.eA = e_org

        fig,ax = plt.subplots()

        ax.plot(e_vec, p_vec_b,label='p')
        ax.set_xlabel(r'$e_{A}$')
        ax.set_ylabel('p*')
        ax.set_title('price')

        plt.show()

    def plot_p_across_e(self):

        par = self.par
        sol = self.sol

        e_vec = np.linspace(3,20,25)
        p_vec_b = np.empty_like(e_vec)
        p_vec_s = np.empty_like(e_vec)

        e_org = par.eA
        for i, e in enumerate(e_vec):

            # solve for e varying for buyer
            par.eA = e
            self.solve(print_sol=False)
            p_vec_b[i] = sol.p

        # reset 
        par.eA = e_org

        e_org = par.eB
        for i, e in enumerate(e_vec):

            # solve for e varying for seller
            par.eB = e
            self.solve(print_sol=False)
            p_vec_s[i] = sol.p

        # reset
        par.eB = e_org

        fig,ax = plt.subplots(ncols=2, sharey=True)
        
        ax[0].plot(e_vec, p_vec_b)
        ax[0].set_xlabel(r'$e_{A}$')
        ax[0].set_ylabel('p*')

        ax[1].plot(e_vec, p_vec_s)
        ax[1].set_xlabel(r'$e_{B}$')
        
        plt.show()