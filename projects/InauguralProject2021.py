import numpy as np
from scipy import optimize

def u_func(c, h, mp):
    """ Calculates utility of chosen (consumption, housing) bundle.

    Args:

    c (float): consumption
    h (float): housing 
    mp (dict): model parameters. 

    Returns:

    (float): utility of bundle
    """
    return (c**(1-mp['phi']))*(h**mp['phi'])

def tau(h, mp, p = 1):
    """ Calculates total housing taxes when choosing h

    Args:

        h (float): housing
        mp (dict): model parameters. 
        p (float): price index of housing

    Returns:

        (float): total taxes paid for a house of quality h

    """

    # Calculate assessment of home. Equation (2).
    p_tilde = p*h*mp['epsilon']
    return mp['tau_g']*p_tilde + mp['tau_p']*(max(p_tilde - mp['p_bar'], 0))

def user_cost(h, mp, p=1):
    """ Get total usercosts of housing, taxes and mortgage payments. Equation (4)

    Args:

        h (float): housing 
        mp (dict): model parameters. 
        p (float): price index of housing

    Returns:

        (float): total user costs of housing. 

    """
    taxes = tau(h, mp, p)
    interest = mp['r']*h*p

    return interest + taxes

def choose_c(h, m, mp, p=1):
    """ Implicit choice of consumption given housing choice. Derived from Equation (3).

    Args:
        h (float): housing 
        m (float): cash-on-hand
        mp (dict): model parameters. 
        p (float): price index of housing

    Returns:

        (float) : consumption given choice of housing and budget constraint. 

    """

    return m - user_cost(h, mp, p)

def value_of_choice(h, m, mp, p=1):
    """ Criterion function for optimizer.

    Args:
        
        h (float): housing 
        m (float): cash-on-hand
        mp (dict): model parameters. 
        p (float): price index of housing

    Returns:

        (float): negative of utility function at (c,h) consumption bundle and cash-on-hand.

    """

    c = choose_c(h, m, mp, p)
    return -u_func(c, h, mp)

def solve_housing(m, mp, print_sol=True, p=1):
    """ Solve the consumers problem given cash-on-hand and model parameters

    Args:
        mp (dict): model parameters. 
        m (float): cash-on-hand
        print_sol (bool): print solution to console
        p (float): price index of housing

    Returns:

        c (float): optimal consumption
        h (float): optimal housing
        u (float): utility at solution

    """

    # Call optimizer  
    sol = optimize.minimize_scalar(value_of_choice, bounds=None,
                                args=(m, mp, p))

    if print_sol:
        print_solution(sol, m, mp, p)

    # Unpack solution
    h = sol.x
    c = choose_c(h, m, mp, p)
    u = u_func(c, h, mp)
    return c, h, u


def tax_revenues(mp, ms, p=1):
    """ Calculates the tax revenue associated with each consumer in the population and its optimal housing

    Args:
        mp (dict): model parameters. 
        ms (np.array): cash-on-hand
        p (float): price index of housing

    Returns:

        (float): distribution of collected housing tax revenue

    """

    h_star = np.empty((len(ms),))
    tax_revenue = np.empty((len(ms),))

    for i,m in enumerate(ms):
        c, h, u = solve_housing(m, mp, print_sol=False, p=p)
        h_star[i] = h
        tax_revenue[i] = tau(h, mp)

    return tax_revenue, h_star

def print_solution(sol, m, mp, p=1):
    """ Print solution of consumer problem

    Args:
        sol (OptimizeResult): solution object from scipy.optimize
        m (float): cash-on-hand
        mp (dict): model parameters.

    Returns:

    """

    h = sol.x
    c = choose_c(h, m, mp, p)
    u = u_func(c, h, mp)

    # Print
    print(f'c          = {c:6.3f}')
    print(f'h          = {h:6.3f}')
    print(f'user_costs = {user_cost(h, mp, p):6.3f}')
    print(f'u          = {u:6.3f}')
    print(f'm - user_costs - c = {m - user_cost(h, mp, p) - c:.4f}')
