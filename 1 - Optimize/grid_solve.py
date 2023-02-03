# All modules used within a module must be imported locally
import numpy as np

# You need to respecify the u_func, because the module does not share scope with the notebook. 
# That is, the module functions cannot see that u_func was defined in the notebook when find_best_choice is called
def u_func(x1,x2,alpha=0.50):
    return x1**alpha * x2**(1-alpha)

def find_best_choice(alpha,I,p1,p2,N1,N2,do_print=True):
    
    # a. allocate numpy arrays
    shape_tuple = (N1,N2)
    x1_values = np.empty(shape_tuple)
    x2_values = np.empty(shape_tuple)
    u_values = np.empty(shape_tuple)
    
    # b. start from guess of x1=x2=0
    x1_best = 0
    x2_best = 0
    u_best = u_func(0,0,alpha=alpha)
    
    # c. loop through all possibilities
    for i in range(N1):
        for j in range(N2):
            
            # i. x1 and x2 (chained assignment)
            x1_values[i,j] = x1 = (i/(N1-1))*I/p1
            x2_values[i,j] = x2 = (j/(N2-1))*I/p2
            
            # ii. utility
            if p1*x1 + p2*x2 <= I: # u(x1,x2) if expenditures <= income 
                u_values[i,j] = u_func(x1,x2,alpha=alpha)
            else: # u(0,0) if expenditures > income, not allowed
                u_values[i,j] = u_func(0,0,alpha=alpha)
            
            # iii. check if best sofar
            if u_values[i,j] > u_best:
                x1_best = x1_values[i,j]
                x2_best = x2_values[i,j] 
                u_best = u_values[i,j]
    
    # d. print
    if do_print:
        print_solution(x1_best,x2_best,u_best,I,p1,p2)

    return x1_best,x2_best,u_best,x1_values,x2_values,u_values

# function for printing the solution
def print_solution(x1,x2,u,I,p1,p2):
    print(f'x1 = {x1:.4f}')
    print(f'x2 = {x2:.4f}')
    print(f'u  = {u:.4f}')
    print(f'I-p1*x1-p2*x2 = {I-p1*x1-p2*x2:.8f}')
    print(f'x1*p1/I = {x1*p1/I:.4f}')


def find_best_choice_monotone(alpha,I,p1,p2,N,do_print=True):
    
    # a. allocate numpy arrays
    shape_tuple = (N)
    x1_values = np.empty(shape_tuple)
    x2_values = np.empty(shape_tuple)
    u_values = np.empty(shape_tuple)
    
    # b. start from guess of x1=x2=0
    x1_best = 0
    x2_best = 0
    u_best = u_func(0,0,alpha)
    
    # c. loop through all possibilities
    for i in range(N):
        
        # i. x1
        x1_values[i] = x1 = i/(N-1)*I/p1
        
        # ii. implied x2 by budget constraint
        x2_values[i] = x2 = (I-p1*x1)/p2
            
        # iii. utility    
        u_values[i] = u_func(x1,x2,alpha)
        
        if u_values[i] >= u_best:    
            x1_best = x1_values[i]
            x2_best = x2_values[i] 
            u_best = u_values[i]
            
    # d. print
    if do_print:
        print_solution(x1_best,x2_best,u_best,I,p1,p2)   

    return x1_best,x2_best,u_best,x1_values,x2_values,u_values