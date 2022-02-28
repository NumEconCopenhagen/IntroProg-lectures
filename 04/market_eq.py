# All modules used within a module must be imported locally
import numpy as np

def demand_good_1_func(alpha,p1,p2,k):
    I = k*p1+p2
    return alpha*I/p1

def demand_good_2_func(alpha,p1,p2,k):
    I = k*p1+p2
    return (1-alpha)*I/p2

def excess_demand_good_1_func(alphas,p1,p2,k):
    
    # a. demand
    demand = np.sum(demand_good_1_func(alphas,p1,p2,k))
    
    # b. supply
    supply = k*alphas.size
    
    # c. excess demand
    excess_demand = demand-supply
    
    return excess_demand

def excess_demand_good_2_func(alphas,p1,p2,k):
    
    # a. demand
    demand = np.sum(demand_good_2_func(alphas,p1,p2,k))
    
    # b. supply
    supply = alphas.size
    
    # c. excess demand
    excess_demand = demand-supply
    
    return excess_demand

def find_equilibrium(alphas,p1_guess,p2,k,kappa=0.5,eps=1e-8,maxiter=500):
    
    t = 0
    p1 = p1_guess
    
    # using a while loop as we don't know number of iterations a priori
    while True:

        # a. step 1: excess demand
        Z1 = excess_demand_good_1_func(alphas,p1,p2,k)
        
        # b: step 2: stop?
        if  np.abs(Z1) < eps or t >= maxiter:
            print(f'{t:3d}: p1 = {p1:12.8f} -> excess demand -> {Z1:14.8f}')
            break    
    
        # c. step 3: update p1
        p1 = p1 + kappa*Z1/alphas.size
            
        # d. step 4: print only every 25th iteration using the modulus operator 
        if t < 5 or t%25 == 0:
            print(f'{t:3d}: p1 = {p1:12.8f} -> excess demand -> {Z1:14.8f}')
        elif t == 5:
            print('   ...')
            
        t += 1    

    return p1