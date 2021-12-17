from numba import njit

@njit
def y_eq_func(ylag,pilag,v,s,slag,alpha,h,b,phi,gamma):
    """ equilibrium value for output

    Args:

        ylag (float): lagged output
        pilag (float): lagged inflation
        v (float): demand disturbance
        s (float): supply disturbance
        slag (float): lagged supply disturbance
        alpha (float): sensitivity of demand to real interest rate
        h (float): coefficient on inflation in Taylor rule
        b (float): coefficient on output in Taylor rule
        phi (float): degree of stickiness in inflation expectations
        gamma (float): effect of output on inflation in SRAS

    Returns:

        (float):  equilibrium value for output

    """

    return 1/(alpha*b+alpha*gamma*h+1)*(-pilag*alpha*h+alpha*gamma*h*phi*ylag+alpha*h*phi*slag-alpha*h*s+v)

@njit
def pi_eq_func(ylag,pilag,v,s,slag,alpha,h,b,phi,gamma):
    """ equilibrium value for inflation

    Args:

        ylag (float): lagged output
        pilag (float): lagged inflation
        v (float): demand disturbance
        s (float): supply disturbance
        slag (float): lagged supply disturbance
        alpha (float): sensitivity of demand to real interest rate
        h (float): coefficient on inflation in Taylor rule
        b (float): coefficient on output in Taylor rule
        phi (float): degree of sticiness in inflation expectations
        gamma (float): effect of output on inflation in SRAS

    Returns:

        (float):  equilibrium value for inflation

    """

    return 1/(alpha*h)*(v-1/(alpha*b+alpha*gamma*h+1)*(alpha*b+1)*(-pilag*alpha*h+alpha*gamma*h*phi*ylag+alpha*h*phi*slag-alpha*h*s+v))

@njit
def _simulate(y,pi,v,s,x_raw,c_raw,alpha,h,b,phi,gamma,delta,omega,sigma_x,sigma_c,T):
    """ equilibrium value for output

    Args:

        y (ndarray): output
        pi (ndarray): inflation
        v (ndarray): demand disturbance
        s (ndarray): supply disturbance
        x_raw (ndarray): raw demand shock
        c_raw (ndarray): raw supply shock
        alpha (float): sensitivity of demand to real interest rate
        h (float): coefficient on inflation in Taylor rule
        b (float): coefficient on output in Taylor rule
        phi (float): degree of sticiness in inflation expectations
        gamma (float): effect of output on inflation in SRAS
        delta (float): persistence of demand shocks
        omega (float): persistence of supply shocks
        sigma_x (float): std. of demand shocks
        sigma_x (float): std. of supply shocks

    """

    for t in range(1,T):
        
        # a. shocks
        v[t] = delta*v[t-1] + sigma_x*x_raw[t]
        s[t] = omega*s[t-1] + sigma_c*c_raw[t]
        
        # b. output
        y[t] = y_eq_func(y[t-1],pi[t-1],v[t],s[t],s[t-1],alpha,h,b,phi,gamma)
        
        # c. inflation
        pi[t] = pi_eq_func(y[t-1],pi[t-1],v[t],s[t],s[t-1],alpha,h,b,phi,gamma)

def simulate(par,sim,T):
    _simulate(sim['y'],sim['pi'],sim['v'],sim['s'],sim['x_raw'],sim['c_raw'],
              par['alpha'],par['h'],par['b'],par['phi'],par['gamma'],
              par['delta'],par['omega'],par['sigma_x'],par['sigma_c'],T)        