from types import SimpleNamespace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
plt.rcParams.update({"axes.grid":True,"grid.color":"black","grid.alpha":"0.25","grid.linestyle":"--"})


class CareerChoiceModel:

    def __init__(self,**kwargs):
        '''
        Initialize model with parameters
        '''
            
        self.par = par = SimpleNamespace()
        self.sim = sim = SimpleNamespace()

        

        for key,value in kwargs.items():
            setattr(par,key,value)

        self.allocate()


    def allocate(self):
        '''
        Allocate size of arrays for simulation
        '''

        par = self.par
        sim = self.sim

        sim.eps = np.empty((par.N,par.J,par.K))
        sim.choice = np.empty((par.N,par.K),dtype=int)
        sim.choice2 = np.empty((par.N,par.K),dtype=int)
        
        
        sim.uf = np.empty((par.N,par.J,par.K))
        sim.u_ante = np.empty((par.N,par.K))
        sim.u_chosen = np.empty((par.N,par.K))
        sim.u_ante2 = np.empty((par.N,par.K))
        sim.u_chosen2 = np.empty((par.N,par.K))

        sim.u_post = np.empty((par.N,par.J,par.K))

        sim.uf2 = np.empty((par.N,par.J,par.K))


    def simulate(self):
        '''
        Simulate shocks and calculate utility of all choices
        '''

        # Draw epsilon:
        par = self.par
        sim = self.sim

        sim.eps[:] = np.random.normal(0,scale=par.sigma, size=sim.eps.shape)

        # Calculate utility
        sim.u_post[:,:,:] = par.v.reshape((1,par.J,1)) + sim.eps

        # draw friends
        for i in range(par.N):
            sim.uf[i,:,:] = par.v.reshape((1,par.J,1)) + np.random.normal(0,scale=par.sigma, size=(par.J,par.F[i],par.K)).mean(axis=1)
        

        
    def utility(self):
        '''
        Determine career choice and store chosen utility
        '''
        par = self.par
        sim = self.sim

        
        sim.choice[:] = np.argmax(sim.uf, axis=1)


        sim.u_ante[:,:] = np.take_along_axis(sim.uf, np.expand_dims(sim.choice, axis=1), axis=1).reshape((par.N,par.K)) # Expected utility from choice
        sim.u_chosen[:] = np.take_along_axis(sim.u_post, np.expand_dims(sim.choice, axis=1), axis=1).reshape((par.N,par.K)) #  Realized utility from choice

    
    def career_change(self):
        '''
        Calculate utility of changing jobs and whether it is optimal
        '''

        sim = self.sim
        par = self.par

        # Expected utility for each job with cost of changing
        sim.uf2[:] = sim.uf - par.c

        # Replace with known utility for the jobs they have chosen
        for j in range(par.J):
            sim.uf2[:,j,:][sim.choice==j] = sim.u_chosen[sim.choice==j]

        
        sim.choice2[:] = np.argmax(sim.uf2, axis=1)
        

        sim.u_ante2[:,:] = np.take_along_axis(sim.uf2, np.expand_dims(sim.choice2, axis=1), axis=1).reshape((par.N,par.K)) # Expected utility from choice
        sim.u_chosen2[:] = np.take_along_axis(sim.u_post, np.expand_dims(sim.choice2, axis=1), axis=1).reshape((par.N,par.K)) - par.c*(sim.choice2 != sim.choice) #  Realized utility from choice



    def plot_post_utility(self):
        par = self.par
        sim = self.sim

        fig, ax = plt.subplots()
        
        
        mean_u_post = sim.u_post.mean(axis=(0,2))


        ax.scatter(range(par.J), mean_u_post,label='Average realized utility')
        ax.scatter(range(par.J), par.v, label='Expected utility',s=8)          


        ax.set_title('Utility of jobs')
        ax.set_xlabel('Job')
        ax.set_ylabel('Utility')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend()


    def plot_utility(self):
        par = self.par
        sim = self.sim


        # Mean utility vs number of friends
        mean_u_ante = sim.u_ante.mean(axis=1)
        mean_u_chosen = sim.u_chosen.mean(axis=1)
        

        fig, ax = plt.subplots()

        
        x_vec = np.arange(1,par.N+1)
        ax.scatter(x_vec, mean_u_ante, label='Expected utility',color='purple')
        ax.scatter(x_vec, mean_u_chosen, label='Realized utility',color='red')
        
        
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title('Mean Utility vs number of friends')
        ax.set_xlabel('Number of Friends')
        ax.set_ylabel('Mean Utility')

        ax.legend()


        # Share choosing each job vs number of friends
        fig, ax = plt.subplots()

        
        
        bottom = np.zeros(par.N)
        for j in range(par.J):
            share = (sim.choice==j).mean(axis=1)*100
            ax.bar(x_vec, share, bottom=bottom, label=f'Share choosing job {j+1}')
            bottom += share
        
        
        ax.set_title('Share choosing each job vs number of friends')
        ax.set_xlabel('Number of Friends')
        ax.set_ylabel('Share, %')

        ax.legend()


    def plot_utility_change(self):
        sim = self.sim
        par = self.par

        
        # Mean utility vs number of friends
        mean_u_ante = sim.u_ante2.mean(axis=1)
        mean_u_chosen2 = sim.u_chosen2.mean(axis=1)

        mean_u_chosen = sim.u_chosen.mean(axis=1)

        fig, ax = plt.subplots()
        
        x_vec = np.arange(1,par.N+1)
        ax.scatter(x_vec, mean_u_ante, label='Expected utility',color='purple')
        ax.scatter(x_vec, mean_u_chosen2, label='Realized utility',color='red')
        ax.scatter(x_vec, mean_u_chosen, label='Chosen utility in earlier period',color='pink')
        
        
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title('Mean Utility vs number of friends')
        ax.set_xlabel('Number of Friends')
        ax.set_ylabel('Mean Utility')

        ax.legend()
        
        
        # Share that change jobs conditional on job chosen in the first period

        fig, ax = plt.subplots()
        
        for j in range(par.J):
            share = np.average(sim.choice2 != sim.choice,axis=1, weights= (sim.choice==j))  *100
            ax.scatter(x_vec, share, label=f'Share initially in job {j+1} that change jobs')
        ax.set_title('Share that change jobs vs number of friends')
        ax.set_xlabel('Number of Friends')
        ax.set_ylabel('Share, %')
        ax.legend()
        


        fig, ax = plt.subplots()
        
        
        share = np.average(sim.choice2 != sim.choice,axis=1)  *100
        ax.scatter(x_vec, share)
        ax.set_title('Share that change jobs vs number of friends')
        ax.set_xlabel('Number of Friends')
        ax.set_ylabel('Share, %')
