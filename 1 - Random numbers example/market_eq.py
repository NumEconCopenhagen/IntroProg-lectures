# All modules used within a module must be imported locally
import numpy as np

class MarketEq():

    def __init__(self,**kwargs):
        '''
        Initialize the model with parameters
        '''

        self.N = 1000 # number of agents
        self.k = 2 # relative endowment of good 1
        self.mu_low = 0.1 # lower bound on alpha
        self.mu_high = 0.9 # upper bound on alpha
        self.kappa = 0.1 # Adjustment factor for solving
        self.eps = 1e-8 # Tolerance parameter for solving
        self.maxiter=500 # Max iterations when solving


        # Update values according to kwargs
        for key, value in kwargs.items():
            setattr(self,key,value) # like self.key = value

        # Simulate alphas according to parameters
        self.simulate_agents()

    def simulate_agents(self):
        '''
        Simulate alphas for all agents
        '''
        self.alphas = np.random.uniform(low=self.mu_low, high=self.mu_high, size=self.N)

    def demand_good_1_func(self,p1,p2):
        I = self.k*p1+p2
        return self.alphas*I/p1

    def demand_good_2_func(self,p1,p2):
        I = self.k*p1+p2
        return (1-self.alphas)*I/p2

    def excess_demand_good_1_func(self,p1,p2):
        
        # a. demand
        demand = np.sum(self.demand_good_1_func(p1,p2))
        
        # b. supply
        supply = self.k*self.N
        
        # c. excess demand
        excess_demand = demand-supply
        
        return excess_demand

    def excess_demand_good_2_func(self,p1,p2):
        
        # a. demand
        demand = np.sum(self.demand_good_2_func(p1,p2))
        
        # b. supply
        supply = self.N
        
        # c. excess demand
        excess_demand = demand-supply
        
        return excess_demand

    def find_equilibrium(self,p1_guess,p2):
        
        t = 0
        p1 = p1_guess
        
        # using a while loop as we don't know number of iterations a priori
        while True:

            # a. step 1: excess demand
            Z1 = self.excess_demand_good_1_func(p1,p2)
            
            # b: step 2: stop?
            if  np.abs(Z1) < self.eps or t >= self.maxiter:
                print(f'{t:3d}: p1 = {p1:12.8f} -> excess demand -> {Z1:14.8f}')
                break    
            
            # c. Print the first 5 and every 25th iteration using the modulus operator 
            if t < 5 or t%25 == 0:
                print(f'{t:3d}: p1 = {p1:12.8f} -> excess demand -> {Z1:14.8f}')
            elif t == 5:
                print('   ...')
            
            # d. step 3: update p1
            p1 = p1 + self.kappa*Z1/self.N
            
            # e. step 4: update counter and return to step 1
            t += 1    


        # Check if solution is found 
        if np.abs(Z1) < self.eps:
            # Store equilibrium prices
            self.p1_star = p1 
            self.p2_star = p2

            # Store equilibrium excess demand 
            self.Z1 = Z1
            self.Z2 = self.excess_demand_good_2_func(self.p1_star, self.p2_star)

            # Make sure that Walras' law is satisfied
            if not np.abs(self.Z2)<self.eps:
                print('The market for good 2 was not cleared')
                print(f'Z2 = {self.Z2}')

        else:
            print('Solution was not found')


    def print_solution(self):

        text = 'Solution to market equilibrium:\n'
        text += f'p1 = {self.p1_star:5.3f}\np2 = {self.p2_star:5.3f}\n\n'

        text += 'Excess demands are:\n'
        text += f'Z1 = {self.Z1}\n'
        text += f'Z2 = {self.Z2}'
        print(text)