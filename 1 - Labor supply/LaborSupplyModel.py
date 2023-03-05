import numpy as np
from scipy import optimize

def implied_tax(l,w,tau0,tau1,kappa):
    """ calculate implied tax of labor supply choice
    
    Args:
    
        l (float): labor supply
        w (float): wage
        tau0 (float): standard labor tax
        tau1 (float): top bracket labor income tax
        kappa (float): cut-off for the top labor income bracket
        
    Returns:
    
        (float): total tax bill
    
    """
    
    return tau0*w*l + tau1*np.fmax(w*l-kappa,0)

def implied_c(l,m,w,tau0,tau1,kappa):
    """ calculate implied optimal consumption of labor supply choice
    
    Args:
    
        l (float): labor supply
        m (float): cash-on-hand
        w (float): wage
        tau0 (float): standard labor tax
        tau1 (float): top bracket labor income tax
        kappa (float): cut-off for the top labor income bracket
        
    Returns:
    
        (float): consumption
    
    """
    
    return m + w*l - implied_tax(l,w,tau0,tau1,kappa)

def utility(c,l,nu,frisch):
    """ utility of consumption and labor supply decision
    
    Args:
    
        c (float): consumption
        l (float): labor supply
        nu (float): disutility of labor supply
        frisch (float): frisch elasticity of labor supply
        
    Returns:
    
        (float): utility
    
    """
    
    return np.log(c) - nu*l**(1+1/frisch)/(1+1/frisch)

def value_of_choice(l,nu,frisch,m,w,tau0,tau1,kappa):
    """ calculate implied utlity of consumption and labor supply choice
    
    Args:
    
        l (float): labor supply
        nu (float): disutility of labor supply
        frisch (float): frisch elasticity of labor supply        
        m (float): cash-on-hand
        w (float): wage
        tau0 (float): standard labor tax
        tau1 (float): top bracket labor income tax
        kappa (float): cut-off for the top labor income bracket
        
    Returns:
    
        (float): utility
        
    """
    
    c = implied_c(l,m,w,tau0,tau1,kappa)
    return utility(c,l,nu,frisch)

def find_optimal_labor_supply(nu,frisch,m,w,tau0,tau1,kappa):
    """ find optimal labor supply choice
    
    Args:
    
        nu (float): disutility of labor supply
        frisch (float): frisch elasticity of labor supply        
        m (float): cash-on-hand
        w (float): wage
        tau0 (float): standard labor tax
        tau1 (float): top bracket labor income tax
        kappa (float): cut-off for the top labor income bracket
        
    Returns:
    
        (float): optimal labor supply
        
    """
    
    obj = lambda l: -value_of_choice(l,nu,frisch,m,w,tau0,tau1,kappa)
    res = optimize.minimize_scalar(obj,bounds=(0,1),method='bounded')

    return res.x