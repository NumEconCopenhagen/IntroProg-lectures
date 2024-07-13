

import numpy as np
from types import SimpleNamespace
import matplotlib.pyplot as plt
from scipy import optimize
import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})



class CO2Model():
    def __init__(self,par):
        '''
        Initialize the model 
        '''
        
        self.par = par
        par.p_vec = np.linspace(0.1,2.0,10)

        # Set parameters for excess_step_solve:
        par.tol = 1e-8
        par.maxiter = 1000
        par.step = 1.

        self.sol = SimpleNamespace()


    def firms(self,p1,p2):

        par = self.par
        sol = self.sol


        w = 1.0 # normalization

        sol.l1 = (p1*par.A*par.gamma/w)**(1/(1-par.gamma))
        sol.y1 = par.A*sol.l1**par.gamma

        sol.l2 = (p2*par.A*par.gamma/w)**(1/(1-par.gamma))
        sol.y2 = par.A*sol.l2**par.gamma

        sol.pi1 = (1-par.gamma)/par.gamma*w*(p1*par.A*par.gamma/w)**(1/(1-par.gamma))
        sol.pi2 = (1-par.gamma)/par.gamma*w*(p2*par.A*par.gamma/w)**(1/(1-par.gamma))
        


    def consumption(self,l,p1,p2):
        par = self.par
        sol = self.sol

        w = 1.0 # normalization

        m = w*l+par.T+sol.pi1+sol.pi2

        sol.c1 = par.alpha*m/p1
        sol.c2 = (1-par.alpha)*m/(p2+par.tau)

        

  
    def households(self,p1,p2):
        par = self.par
        sol = self.sol

        def value_of_choice(l):
            
            self.consumption(l[0],p1,p2)
            u_c = np.log(sol.c1**par.alpha*sol.c2**(1-par.alpha))
            u_l = par.nu*l[0]**(1+par.epsilon)/(1+par.epsilon)

            return -(u_c-u_l)

        res = optimize.minimize(value_of_choice,0.5,bounds=((0.0,None),),
                                method='nelder-mead')
        
        sol.l = res.x[0]
        sol.U = -res.fun
        

    def market_clearing(self,p):
        par = self.par
        sol = self.sol

        p1,p2 = p
        self.firms(p1,p2)
        # Setting T here, is a somewhat subtle point. 
        # It is important to set T before solving the household problem as it affects the household solution.
        # It can be set here because equlibrium implies that c2=y2 so only the firm problem is needed to find T.
        par.T = par.tau*sol.y2

        self.households(p1,p2)

        return (sol.y1-sol.c1),(sol.y2-sol.c2),(sol.l-sol.l1-sol.l2)



    def solve_grid(self):
        '''
        Solve for approximate equilibrium prices over a grid
        
        '''

        par = self.par 
        sol = self.sol


        check = np.inf
        p1_ast = None
        p2_ast = None

        for p1 in par.p_vec:
            for p2 in par.p_vec:

                errors = self.market_clearing((p1,p2))
                check_now = errors[0]**2+errors[1]**2
                
                if check_now < check:

                    check = check_now
                    p1_ast = p1
                    p2_ast = p2



        print(f'{p1_ast = :.2f}, {p2_ast = :.2f}')

        sol.p1_approx = p1_ast
        sol.p2_approx = p2_ast


    def solve(self,guess='gridsol',method='minimize',do_print=False):
        '''
        Solve for equilibrium prices
        '''

        par = self.par
        sol = self.sol

        def obj(p):
            errors = self.market_clearing(p)
            return errors[0],errors[1]
        
        if guess == 'gridsol':
            guess = (sol.p1_approx,sol.p2_approx)
        elif guess == 'sol':
            guess = (sol.p1,sol.p2)            
        
        if method in ['minimize']:
            res = optimize.root(obj, guess, method='hybr')
            sol.p1 = res.x[0]
            sol.p2 = res.x[1]
        elif method in ['excess_step']:
            sol.p1,sol.p2 = self.excess_step_solve(guess,do_print)
           

        if do_print:
            print(f'{sol.p1 = :.2f}, {sol.p2 = :.2f}')

            errors = self.market_clearing((sol.p1,sol.p2))
            print(f'{errors = }')
        

    def excess_step_solve(self,guess,do_print):
        '''
        Find equilibrium prices using excess demand to determine the step size
        '''
        par = self.par
        p1,p2 = guess
        t = 0 
        
        while True:

            # exess demand
            error1, error2,error3 = self.market_clearing((p1,p2))
            
            # stop? 
            if np.abs(error1) < par.tol and np.abs(error2) < par.tol:
                if do_print:
                    print(f'{t:3d}: p1 = {p1:6.8f}, p2 = {p2:6.8f} -> excess goods demand -> {error1:14.8f} {error2:14.8f}')
                break
            
            # Print 
            if do_print:
                if t < 5 or t%25 == 0:
                    print(f'{t:3d}: p1 = {p1:6.8f}, p2 = {p2:6.8f} -> excess goods demand -> {error2:14.8f} {error2:14.8f}')
                elif t == 5:
                    print('   ...')


            # Update prices
            
            p1 -= par.step*error1
            p2 -= par.step*error2
            

            # No more iterations? 
            if t >= par.maxiter:
                print('Solution not found!')
                break

            t += 1
        
        return p1,p2


    def optimal_gov(self):
        '''
        Find optimal values of tau and T for the government
        '''
        par = self.par
        sol = self.sol

        def value_of_choice_gov(tau,do_print=False):
            par = self.par
            sol = self.sol
            
            par.tau = tau
            
            self.solve(guess='sol')

            self.firms(sol.p1,sol.p2)
            par.T = par.tau*sol.y2
            self.households(sol.p1,sol.p2)

            SWF = sol.U - par.kappa*sol.y2    

            if do_print: print(f'{tau = :.3f}, {sol.U = :.4f}, {sol.y2 = :.4f}, {SWF = :.5f}')

            return -SWF    


        taus = np.linspace(0.0,0.3,20)
        SWFs = np.zeros_like(taus)

        print('Grid search:')
        for i,tau in enumerate(taus):
            SWFs[i] = -value_of_choice_gov(tau,do_print=True)

        i = np.argmax(SWFs)
        tau_min = taus[i-1]
        tau_max = taus[i+1]

        res = optimize.minimize_scalar(value_of_choice_gov,bounds=(tau_min,tau_max),method='bounded')

        print('\nOptimal tau:')
        value_of_choice_gov(res.x,do_print=True)
        T = par.tau*sol.y2
        print(f'Gvining optimal T = {T:.2f}')


        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(taus,SWFs,label='SWF');
        ax.axvline(res.x,ls='--',color='k',label=r'Optimal $\tau$')

        ax.set_xlabel(r'$\tau$')
        ax.legend()
        fig.tight_layout()