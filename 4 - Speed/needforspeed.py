import numpy as np

def create_list(n):
    return [i*i for i in range(n)]

def create_generator(n):
    return (i*i for i in range(n))

def test_memory(n):

    # list vs. generators
    squares = create_list(n)
    squares_generator = create_generator(n)

    # numpy arrays of different types
    A = np.ones((1000,1000))
    B = np.ones((1000,1000),dtype=np.double)
    C = np.ones((1000,1000),dtype=np.single)
    D = np.ones((1000,1000),dtype=np.int64)
    E = np.ones((1000,1000),dtype=np.int32)
    F = np.ones((1000,1000),dtype=np.bool)