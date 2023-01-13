from math import isnan
from math import fabs
import numpy as np


def safe_NR(f, fprime, a, b, x0=np.nan, maxit=20, tol=1e-14, do_print=False):
    """ 
    Root of f(x) obtained by switching between Newton-Raphson and bisection.
    Args:
        f (function): function to find root of
        fprime (function): derivative of f
        a (float): lower bound
        b (float): upper bound
        x0 (float): initial guess of root
        maxit (int): maximum number of iterations
        tol (float): tolerance
    """  
    """
    Numerical root of f(x) = 0 by Newton-Raphson method
    """

    if a > b :
        raise Exception('Lower bound exceeds upper bound')

    if f(a)*f(b) > 0:
        raise Exception(
            'There is either no root in interval or multiple roots. sign(f(a)) == sign(f(b))')

    if (fabs(f(a)) < tol):
        return a 

    if (fabs(f(b)) < tol):
        return b 

    # Initial values 
    if isnan(x0):
        x = (b+a)/2
    else:
        x = x0
  
    fx = f(x)
    dfx = fprime(x)
    j = 0 

    NRstep = lambda x,fx,dfx : x - fx/dfx
    Bstep = lambda a,b: a + 0.5*(b-a)  

    while j < maxit:
        cond1 = not (a < NRstep(x,fx,dfx) < b)
        if cond1:
            # Go Bisection if Newton-Raphson goes off bounds or goes too slow 
            xn = Bstep(a,b)
            steptype = 'bisection'
        else:
            # Default to Newton-Raphson
            xn = NRstep(x,fx,dfx)
            steptype = 'newton'

        if do_print:
            print(j, ' a:', a, 'b:', b, 'xn:', xn, 'fxn:', f(xn), ' steptype:', steptype)

        x = xn 
        fx = f(x)
        dfx = fprime(x)

        if fabs(fx) < tol: 
            break

        # Update bounds. Needed for both bisection and Newton-Rapshon  
        if fx*f(a) < 0.0:
            b = x 
        else:
            a = x 
        
        j = j+1

    return x,j 
                