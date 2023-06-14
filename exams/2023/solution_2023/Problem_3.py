import time
from types import SimpleNamespace
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})
plt.rcParams.update({'font.size': 14})

class GlobalOptimizerClass:

    def __init__(self,f,bounds,K,K_ubar,tol=1e-8):
        """ setup """

        self.f = f

        # step 1
        self.bounds = bounds

        # step 2
        self.K = K
        self.K_ubar = K_ubar
        self.tol = tol

    def run(self,do_print=True):
        """ run the algorithm """
        
        xs = np.zeros((self.K,2))
        fs = np.zeros(self.K)

        # step 3
        for k in range(self.K):

            # step A
            xk = np.empty(2)
            xk[0] = np.random.uniform(self.bounds[0,0],self.bounds[0,1])
            xk[1] = np.random.uniform(self.bounds[1,0],self.bounds[1,1])

            # step D
            if not k < self.K_ubar:

                chi_k = 0.5*2/(1+np.exp((k-self.K_ubar)/100))
                xk = chi_k*xk + (1-chi_k)*x_ast

            xs[k] = xk

            # step E
            res = optimize.minimize(self.f,xk,method="BFGS",tol=self.tol)
            fs[k] = res.fun

            # step F
            if k == 0 or res.fun < f_ast:

                if do_print: print(f'k = {k:4d}, f_ast = {res.fun:12.8f}, x_ast = {res.x}')
                    
                f_ast = res.fun
                x_ast = res.x

            # step G
            if res.fun < self.tol: break

        # step 4
        self.res = SimpleNamespace(x_ast=x_ast,f_ast=f_ast,k=k,xs=xs[:k,:],fs=fs[:k])

    def plot(self):
        """ plot the results"""

        res = self.res

        fig = plt.figure(figsize=(12,4))
        ax = fig.add_subplot(1,2,1)
        ax.set_title('effective initial guess')
        ax.scatter(np.arange(res.xs.shape[0]),res.xs[:,0],label='$x_0$')
        ax.scatter(np.arange(res.xs.shape[0]),res.xs[:,1],label='$x_1$')
        #ax.set_yscale('symlog')
        ax.legend()

        ax = fig.add_subplot(1,2,2)
        ax.set_title('objective function')
        ax.scatter(np.arange(res.xs.shape[0]),res.xs[:,1],label='$x_1$')
        #ax.set_yscale('symlog');        

    def test(self,reps=10):
        """ run the algorithm multiple times """

        times = np.empty(reps)
        iterations = np.empty(reps)
        for i in range(reps):

            t0 = time.time()
            self.run(do_print=False)
            res = self.res
            times[i] = time.time()-t0
            iterations[i] = res.k
            print(f'Elapsed time: {times[i]:5.2f} seconds. Iterations: {res.k:4d}. f_ast = {res.f_ast:5.2e}')        


        print(f'Average time:       {np.mean(times):5.2f} seconds')
        print(f'Average iterations: {np.mean(iterations):5.2f}')