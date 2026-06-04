

# standard numerics and linear algebra libraries
import numpy as np  
import numpy.linalg as LA
import scipy.linalg as sciLA

# sparse matrix functions
import scipy.sparse as sparse

# for interactive graphics
from ipywidgets import interactive, interact
from ipywidgets import FloatSlider

# avoid typing np.XY all the time
from numpy import (array, pi, cos, sin, ones, size, sqrt, real, mod, append, arange, exp)
from Adva_Quant_Inform.operators import buildSparseGateSingle, initRegister, indToState, buildUf 

def deutschJosza(f, n):
    psi = initRegister(n)
    # apply the Hadamards
    for i in arange(n):
        psi = buildSparseGateSingle(n,i,H) @ psi
    # apply U_f
    psi = buildUf(f, n) @ psi
    # apply the Hadamards again
    for i in arange(n):
        psi = buildSparseGateSingle(n,i,H) @ psi

    # If the probability of having the all zero state is 1, then f is constant.
    # Since the state of all zero is represented in the computational basis by 1 in the first entry 
    # and then all zeros, one can just check np.isclose(np.abs(psi[0])**2, 1), 
    # namely the probability equal to one for the first element of the vector psi[0].
    # Even simpler, f is constant iff psi[0]=\pm 1 the function is constant (we have a binary choice).
    if psi[0] == 0:
        print('The function is balanced.')
    else:
        print('The function is constant.')

    # checking
    ratio = np.sum([f(indToState(n,x)) for x in range(2**n)])/2**n
    print("The ratio of ones to zeros (computed directly) is:", ratio)
    return psi 
