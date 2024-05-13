import numpy as np
import scipy as sp
import sympy as sm

def gauss_jordan(A,no_reduction=False):
    """ Gauss-Jordan elimination

    Convert a matrix to its row (reduced) echelon form.

    Args:

        A (ndarray): input matrix to work on in-place
        no_reduction (bool): stop when non-reduced echelon form is reached

    """

    n,_m = A.shape

    # row echelon form (Gauss elimination)
    for i in range(0, n):

        # a. search for maximum in this column
        maxrow = i + np.argmax(abs(A[i:,i]))

        # b. swap maximum row with current row (column by column)
        temp = A[maxrow,i:].copy()
        A[maxrow,i:] = A[i,i:]
        A[i,i:] = temp

        # b. make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k,i]/A[i,i]
            A[k,i] = 0
            A[k,i+1:] += c*A[i,i+1:]
        
    if no_reduction:
        return

    # reduced row echelon form (Gauss-Jordan elimination)
    for i in range(n-1,-1,-1):
        
        # a. normalize this row
        c = A[i,i]
        A[i,:] /= c

        # b. make all rows above this one 0 in the current column
        for j in range(0,i):
            c = A[j,i]
            A[j,:] -= c*A[i,:]

def gauss_seidel_split(A):
    """ split A matrix in additive lower and upper triangular matrices

    Args:

        A (ndarray): input matrix

    Returns:

        L (ndarray): lower triangular matrix
        U (ndarray): upper triangular matrix (zero diagonal)

    """

    L = np.tril(A)
    U = np.triu(A)
    np.fill_diagonal(U,0)
    return L,U

def solve_with_forward_substitution(L,RHS):
    """ solve matrix equation with forward substitution

    Args:

        L (ndarray): lower triangular matrix
        RHS (ndarray): vector of right-hand-side variables

    Returns:

        x (ndarray): Solution vector

    """

    n = RHS.size
    x = np.zeros(n)
    for i in range(n):
        x[i] = RHS[i]
        for j in range(i):
            x[i] -= L[i,j]*x[j]    
        x[i] /= L[i,i]
    
    return x

def solve_with_backward_substitution(U,RHS):
    """ solve matrix equation with backward substitution

    Args:

        L (ndarray): uppper triangular matrix
        RHS (ndarray): vector of right-hand-side variables

    Returns:

        x (ndarray): Solution vector

    """

    n = RHS.size
    x = np.zeros(n)
    for i in reversed(range(n)):
        x[i] = RHS[i]
        for j in range(i+1,n):
            x[i] -= U[i,j]*x[j]    
        x[i] /= U[i,i]
    
    return x

def gauss_seidel(A,b,x0,max_iter=500,tau=10**(-8),do_print=False):
    """ solve matrix equation with Gauss-Seidel

    Args:

        A (ndarray): LHS matrix
        b (ndarray): RHS vector
        x0 (ndarray): guess on solution
        max_iter (int): maximum number of iterations (optional)
        tau (float): tolerance level
        do_print (bool): indicator for whether to print or not

    Returns:

        x (ndarray): Solution vector

    """

    converged = False

    # a. split
    L,U = gauss_seidel_split(A)
    
    # b. iterate
    x = x0
    i = 0

    if do_print:
        print('  ',x)

    while i < max_iter and not converged:
        
        # i. save previous
        x_prev = x
        
        # ii. compute RHS
        y = b-U@x

        # iii. solve with forward substituion
        x = solve_with_forward_substitution(L,y)
        #x = sp.linalg.solve_triangular(L,y,lower=True) # equivalent, but faster
        
        # iv. check convergence
        max_abs_diff = np.max(np.abs(x-x_prev))
        if max_abs_diff < tau:
            converged = True

        # v. misc
        if do_print:
            print(f'{i:2d}',x)

        i += 1

    return x

def lu_decomposition(A):
    """ compute LU decomposition

    Args:

        A (ndarray): input matrix

    Returns:

        L (ndarray): lower triangular matrix
        U (ndarray): upper triangular matrix

    """
    
    n = len(A)

    # a. create zero matrices for L and U                                                                                                                                                                                                                 
    L = np.zeros((n,n))
    U = np.zeros((n,n))

    # b. set diagonal of L to one
    np.fill_diagonal(L,1)
    
    # c. perform the LU Decomposition                                                                                                                                                                                                                     
    for j in range(n):          

        for i in range(j+1):
            c = U[:,j]@L[i,:]
            U[i][j] = A[i][j] - c

        for i in range(j, n):
            c = U[:j,j]@L[i,:j]
            L[i][j] = (A[i][j] - c) / U[j][j]

    return L,U

def construct_sympy_matrix(positions,name='a'):
    """ construct sympy matrix with non-zero elements in positions
    
    Args:
    
        Positions (list): list of positions in strings, e.g. ['11','31']
    
    Returns:
    
        mat (sympy.matrix): Sympy Matrix
    
    """
    
    # a. dictionary of element with position as key and a_position as value
    entries = {f'{ij}':sm.symbols(f'{name}_{ij}') for ij in positions}

    # b. function for creating element or zero
    add = lambda x: entries[x] if x in entries else 0

    # c. create matrix
    mat_as_list = [[add(f'{1+i}{1+j}') for j in range(3)] for i in range(3)]
    mat = sm.Matrix(mat_as_list)

    return mat

def fill_sympy_matrix(A_sm,A,name='a'):

    n,m = A.shape

    # a. make all substitution
    A_sm_copy = A_sm
    for i in range(n):
        for j in range(m):
            if not A[i,j] == 0:
                A_sm_copy = A_sm_copy.subs(f'{name}_{1+i}{1+j}',A[i,j])
    
    # b. lambdify with no inputs
    f = sm.lambdify((),A_sm_copy)

    # c. return filled matrix
    return f()
                        
                                    
