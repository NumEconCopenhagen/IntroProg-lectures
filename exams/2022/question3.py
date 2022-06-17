import numpy as np

def f_approx(x, f, N, M):
    '''Chebyshev approximation of the function f at x

    Args:
        x (ndarray)     : the points at which to approximate f
        f (callable)    : function to approximate
        N (int)         : degree of approximation
        M (int)         : number of evaluation points of f

    Returns:
        f_hat (ndarray) : approximation to f at x 
    '''

    # Basis function T
    T = lambda i,x: np.cos(i*np.arccos(x))

    # Step 1: find nodes on function interval    
    k = np.arange(1, M+1)
    zk = -np.cos(np.pi*(2*k-1)/(2*M))
    
    # Step 2: evaluate true function at nodes
    yk = f(zk)

    # Step 3: compute coefficients ai for approximation. ai has dimension (N+1)x1
    ai = np.empty((N+1,1))
    for i in range(N+1):
        ai[i] = (yk @ T(i,zk))/np.sum(T(i,zk)**2)

    # Create matrix of x with repeated values in row dimension for convenient multiplication 
    # X has dimension (N+1) x J, where J is number of elements in x
    x = np.reshape(x,(len(x),1))
    X = np.transpose(np.tile(x, (1,N+1)))
    iN = np.arange(N+1, dtype='int')
    iN = iN[np.newaxis,:]
    
    # Step 4.1: compute T(x) at all i=0,..,N (row-dimension) for all elements in x (column-dimension)
    Tx = T(iN.T, X)

    # Step 4.2: compute approximation by summing the product of coefficients ai and node values Ti(x)
    # Note that ai is multiplied with each column of Tx, which represent a separate point to approximate f at.
    f_hat = np.sum(ai*Tx, axis=0)

    return f_hat
