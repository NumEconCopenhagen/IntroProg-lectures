from operator import le
import numpy as np
from scipy import interpolate
from scipy import optimize

def utility(c,mp):
    ''' Utility function of agent 
    
    Args:
        c (float): consumption
        mp (SimpleNamespace): model parameters
    '''
    return c**(1-mp.rho)/(1-mp.rho)

def bequest(c,m,mp):
    ''' Bequest function of agent. Added to utility in last period.

    Args:
        c (float): consumption
        m (float): cash-on-hand
        mp (SimpleNamespace): model parameters
    
    '''
    return mp.nu*(m-c+mp.kappa)**(1-mp.rho)/(1-mp.rho)

def v_last_period(c,m,mp):
    ''' Value function in last period.
    '''
    return utility(c,mp) + bequest(c,m,mp)

def solve_last_period(mp):
    ''' Solving the model in last period by allocating cash-on-hand into consumption and bequests.
    
    Args:
        mp (SimpleNamespace): model parameters
    
    Returns:
        m_grid (ndarray): applied grid over m
        v_func (ndarray): value function at each m in m_grid
        c_func (ndarray): consumption policy function at each m in m_grid
    '''

    # a. Set up containers for output 
    m_grid = np.linspace(mp.m_min,mp.m_max,mp.Nm)
    v_func = np.empty(mp.Nm)
    c_func = np.empty(mp.Nm)

    # b. Solve last period on grid over m
    for i,m in enumerate(m_grid):

        # i. define objective
        obj = lambda x: -v_last_period(x[0],m,mp)

        # ii. call optimizer
        x0 = m/2 # initial value
        result = optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((1e-8,m),))

        # iii. save results
        v_func[i] = -result.fun
        c_func[i] = result.x
        
    return m_grid,v_func,c_func

def v(c,s,m,mp,v_nxt_interp):
    ''' Expected value function

    Args:
        c (float): consumption
        s (int): school choice
        m (float): cash-on-hand
        mp (SimpleNamespace): model parameters
        v_nxt_interp (callable): interpolator over next period value function

    Retuns:
        (float): choice-specific expected value
    '''
    
    # a. expected value
    v_nxt = 0.0
    prbs = np.array([mp.p, 1-mp.p])
    ys = mp.ybar + mp.gamma*s + np.array([mp.Delta,-mp.Delta]) 

    for prb,y in zip(prbs, ys):
        
        # i. cash-on-hand in next period
        a = m-c
        m_nxt = (1+mp.r)*a + y
        
        # ii. values in next period 
        v_nxt_y = v_nxt_interp([m_nxt])[0]
        
        # iii. probability weighted sum of values
        v_nxt += prb*v_nxt_y
    
    # b. total value
    return utility(c,mp) + mp.beta*v_nxt

def solve_single_period(mp,v_nxt_interp,m_grid=None):
    ''' Solve a single period's objective of choosing optimal consumption and schooling level

    Args:
        mp (SimpleNamespace): model parameters
        v_nxt_interp (callable): interpolator over next period value function

    Retuns:
        m_grid (ndarray): applied grid over m
        v_func (ndarray): value function at each m in m_grid
        c_func (ndarray): consumption policy function at each m in m_grid
        s_func (ndarray): school choice function at each m in m_grid
        v_schoice (ndarray): school choice specific value functions. s=0 in first column, s=1 in second.
    '''

    # a. allocate
    if m_grid is None:
        m_grid = np.linspace(mp.m_min,mp.m_max,mp.Nm)
        
    v_func = np.empty(mp.Nm)
    c_func = np.empty(mp.Nm)
    s_func = np.empty(mp.Nm)
    v_schoice = np.empty((mp.Nm,2))

    # b. solve period problem of choosing s and c
    for i,m in enumerate(m_grid):
        v_s = []
        # i. solve for s=0 and s=1 and store results in v_s
        for s in [0,1]:
            m_s = m - s*mp.tau # school specific cash-on-hand, reduced by tuition when s=1
            obj = lambda x: -v(x[0],s,m_s,mp,v_nxt_interp)
            x0 = m_s/2 # initial guess
            v_s.append(optimize.minimize(obj,[x0],method='L-BFGS-B',bounds=((1e-8,m_s),)))
        
        # ii. take the envelope of values for {s=0, s=1}
        s = int(-v_s[1].fun > -v_s[0].fun)

        # iv. save result
        s_func[i] = s
        v_func[i] = -v_s[s].fun
        c_func[i] = v_s[s].x[0]
        v_schoice[i,:] = [-v_s[0].fun, -v_s[1].fun]
     
    return m_grid, v_func, c_func, s_func, v_schoice


def solve(mp):
    ''' Solve 2 period model of schooling choice and consumption-saving decision by backwards induction

    Args:
        mp (SimpleNamespace): model parameters

    Retuns:
        m1_grid (ndarray): applied grid over initial cash-on-hand
        c1_func (ndarray): consumption policy function in first period at each m in m1_grid
        s_func (ndarray): school choice function at each m in m1_grid
        v1_func (ndarray): value function in first period at each m in m1_grid
        m2_grid (ndarray): applied grid over period 2 cash-on-hand
        c2_func (ndarray): consumption policy function in second period at each m in m2_grid
        v1_schoice (ndarray): school choice specific value functions in first period. s=0 in first column, s=1 in second.
    '''
    
    # a. solve period 2
    m2_grid,v2_func,c2_func = solve_last_period(mp)
    
    # b. construct interpolator
    v2_func_interp = interpolate.RegularGridInterpolator([m2_grid], v2_func,
        bounds_error=False,fill_value=None)
    
    # b. solve period 1
    m1_grid,v1_func,c1_func,s_func,v1_schoice = solve_single_period(mp,v2_func_interp)
    
    return m1_grid, c1_func, s_func, v1_func, m2_grid, c2_func,v2_func, v1_schoice



def bisection_search(pmin, pmax, parname, model_fun, delta, mp, eps = 1e-5, do_print=False):
    """ Search for parameter value that changes a model outcome by delta units.
        REQUIREMENTS: pvals must be sorted and the outcome of model_fun must be monotone in pvals.

        Example: the school choice policy function may change from 0 to 1 when changing the tuition paramter.  

    Args:
    
        pmin (float):           Smallest potential parameter candidate.
        pmax (float):           Largest potential parameter candidate.
        parname (str):          Name of parameter which is being searched over.
        model_fun (callable):   A reference to the model for which the change in outcome is wanted. Must return a float.
        delta (float):          Size of change in outcome that model_fun must produce.
        mp (SimpleNamespace):   Model parameters for model_fun().
        eps (float):            Threshold for nearness to point where model outcome changes.
        do_print (bool):        Indicator for printing progress.

        
    Returns:
    
        pmin (float):   The first parameter value BEFORE model outcome changes by delta. 
        max (float):    The first parameter values AFTER model outcome changes by delta. 
        it (int):       Number of bisection iterations.
    
    """

    # a. initialize values
    found = False

    # b. obtain outcome for minmial parameter guess
    setattr(mp, parname, pmin)
    v_first = model_fun(mp) 
    
    # c. test that first and last candidate parameter value yields a difference equal to delta arg
    setattr(mp, parname, pmax)
    v_last = model_fun(mp) 
    assert np.isclose(np.abs(v_last - v_first), delta), "Difference in outcome for first and last element in pvals does not equal delta input"

    # c. bisection loop
    it = 0
    while pmin < (pmax-eps) and not found :

        # i. find midpoint in remaining interval of potential parameter values
        pmid = pmin + (pmax-pmin)/2.0 

        # ii. update outcome based on new midpoint 
        setattr(mp, 'tau', pmid)
        v_mid = model_fun(mp) 
        
        if do_print:
            print(f'{it:2d} current midpoint = {pmid:3.6f}')

        # iii. if outcome at new midpoint is delta apart from first candidate value, the candidate of interest is between first element and midpoint. 
        # Therefore discard everything to the right. Otherwise, discard everything to the left of midpoint 
        if np.isclose(np.abs(v_mid - v_first), delta):
            pmax = pmid 
        else: 
            pmin = pmid
        
        if pmax-pmin < eps:
            found = True

        it += 1 

    # d. result: pmax and pmin are now eps-distance apart and the kink is somewhere in between
    return pmin,pmax,it
