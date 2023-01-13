import numpy as np

def log_likelihood(theta, x, y, mp):
    ''' Computes the sum of log-likelihood contributions for observations given beta
    
    Args:
        
        beta (ndarray): coefficients to independent variables
        x (ndarray): independent variables
        y (ndarray): observed binary choices
        
    Returns:
    
        (float): sum of log-likelihood contributions
  
    '''
    b1 = theta[0]
    b2 = theta[1] + theta[2] * mp.Z  # Drawing individual slopes
    b2 = np.broadcast_to(b2, (mp.N, mp.Ndraws))

    xb = b1 + x*b2 
    z = np.exp(xb)
    y_prb = z / (1 + z)
    
    # Calculate likelihood contributions
    y_dens = y*y_prb + (1-y)*(1-y_prb)
    y_mc = np.mean(y_dens, axis=1)
    ll = np.sum(np.log(y_mc))  
    return ll
