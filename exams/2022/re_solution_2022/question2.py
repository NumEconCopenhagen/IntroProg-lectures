import numpy as np

def utility(w, s):
    ''' Utility from social media platform usage. 
    ''' 
    return w + np.sum(s, axis=0)


def utility_young(w, s):
    ''' Utility from social media platform usage for young users.
    '''
    return w + s[0] - s[1]


def ccp(u):
    ''' Choice probabilities.
    '''
    return np.exp(u) / np.sum(np.exp(u), axis=1, keepdims=True)


def share(prb, mp):
    return prb[0]*mp.rho_y + prb[1]*mp.rho_o



def some_eq(mp, s0):
    ''' Find user distribution on social media platforms by succesive approximations

    Args:
        mp (SimpleNamespace): model parameters
        s0 (ndarray): initial vector of user shares for platforms

    Returns:
        (ndarray): equilibrium share of each user type on each platform
    '''

    # a.1 Collect population shares in vector         
    rho = np.array([[mp.rho_y],[mp.rho_o]])

    # a.2 Use desired specification of utility for the young and the old
    u_old = lambda s : utility(mp.w_o, s)
    if (mp.utility_spec == 'baseline'):
        u_young = lambda s: utility(mp.w_y,s)
    elif (mp.utility_spec == 'negative_externality'):
        u_young = lambda s: utility_young(mp.w_y,s)
    else:
        raise Exception('wrong utility specification')
    
    # a.3 Compute initial utility vector
    s = np.stack((s0,s0)) * rho
    uk_young = u_young(s)
    uk_old = u_old(s)
    u_prev = np.stack((uk_young,uk_old))
    k = 0

    # b.1 Succesive approximation loop
    while k < mp.N:
        s = ccp(u_prev) * rho
        uk_young = u_young(s)
        uk_old = u_old(s)
        u_nxt = np.stack((uk_young, uk_old))

        if np.amax(np.abs(u_prev - u_nxt)) < mp.eps:            
            break
        else:
            k = k+1
            u_prev = u_nxt

    s_total = np.sum(s, axis=0)
    s_age = s / rho
    return s_total, s_age, u_nxt
