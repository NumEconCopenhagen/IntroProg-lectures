from types import SimpleNamespace
import math
import numpy as np
from scipy import optimize

def util(y,theta=-2):
    return (y**(1+theta))/(1+theta)

def premium(q,p,l=0):
    return (1+l)*p*q

def V(q,x,y,p):
    z = y-x+q-premium(q,p)
    return p*util(z) + (1-p)*util(y-premium(q,p))

def V_tilde(pi,p,q,x,y):
    z = y-x+q-pi
    return p*util(z) + (1-p)*util(y-pi)


def opt_q(x, y, p):
    def crit(q): return -V(q, x, y, p)

    q0 = [x/2]

    sol = optimize.minimize_scalar(
        crit, q0, method='bounded', bounds=(1e-5, x))

    return sol

def mc_value(mp, pol):
    x = np.random.beta(mp.a, mp.b, mp.N)
    z = mp.y - (1-pol.gamma)*x - pol.pi
    val = np.mean(util(z))
    return val

