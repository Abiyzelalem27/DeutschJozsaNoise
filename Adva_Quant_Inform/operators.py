


# standard numerics and linear algebra libraries
import numpy as np  
import numpy.linalg as LA
import scipy.linalg as sciLA

# sparse matrix functions
import scipy.sparse as sparse

# for interactive graphics
from ipywidgets import interactive, interact
from ipywidgets import FloatSlider 
from qutip import tensor, qeye, sigmax, sigmay, sigmaz
# avoid typing np.XY all the time
from numpy import (array, pi, cos, sin, ones, size, sqrt, real, mod, append, arange, exp)

# Define single qubit gates
X = array([[0,1],
           [1,0]])
Y = array([[0,-1j],
           [1j,0]])
Z = array([[1,0],
           [0,-1]])
H = array([[1,1],
           [1,-1]])/sqrt(2)
S = array([[1,0],
           [0,-1j]])
T = array([[1,0],
           [0,exp(1j*pi/4)]])

P0 = np.array([[1,0],
            [0,0]])
P1 = np.array([[0,0],
            [0,1]]) 

I = np.identity(2)

# helper function for generating basis vectors
def basisvec(n, k):
    v = np.zeros(2**n)
    v[k] = 1
    return v 
    
def rotation(ax,theta):
    return sciLA.expm(-1j * theta/2 * (ax[0]*X + ax[1]*Y + ax[2]*Z)) 
    
def buildSparseGateSingle(n, i, gate):
    sgate = sparse.csr_matrix(gate)
    return sparse.kron(sparse.kron(sparse.identity(2**i), sgate), sparse.identity(2**(n-i-1)))

def buildSparseCNOT(n, ic, it):
    P0ic = buildSparseGateSingle(n, ic, P0)
    P1ic = buildSparseGateSingle(n, ic, P1)
    Xit  = buildSparseGateSingle(n, it, X)
    return P0ic + P1ic @ Xit

def initRegister(n):
    return basisvec(n,0)

def indToState(n, k):
    num = bin(k)[2:].zfill(n)
    return array([int(x) for x in str(num)])

def stateToInd(state):
    return int("".join(str(x) for x in state),2)

def buildUf(f, n):
    return sparse.diags([(-1)**f(indToState(n,x)) for x in range(2**n)]) 

def systemSizeFromState(psi):
    return int(np.log2(len(psi)))

# The following function picks a random vector out of the possible outcomes 
def doMeasurement(psi):
    n = systemSizeFromState(psi)
    pvec = np.abs(psi)**2
    thresholds = np.cumsum(pvec)
    r = np.random.rand()
    indOutcome = np.sum(thresholds < r)
    return indToState(n, indOutcome)


# The previous task can be also performed differently.
# In a more abstract language, one might do:

# ind = 0
# cumprob = 0
# r = rand(0,1)
# while cumprob < r
#    cumprob += P(ind)
#    ind++

# where the index of interest when one stops is ind-1. 
# P(ind) is the probability associated to the outcome identified by ind.

def varianceSigmaZ(x):
    return 1 - (2*x**2-1)**2  
# Define the function to maximize (we negate it to use minimize_scalar)
def NegVarianceSigmaZ(x):
    return -varianceSigmaZ(x)

# pre-allocate operators
si = qeye(2) # identity
sx = sigmax()
sy = sigmay()
sz = sigmaz()

# This is basically the same as the function
# buildSparseGateSingle(n, i, gate) from exercise sheet 3 does.
def sx_list(N):
    sx_list_flag = []
    for n in range(N):
        op_list = []
        for m in range(N):
            op_list.append(si)
        op_list[n] = sx
        sx_list_flag.append(tensor(op_list))
    return(sx_list_flag)

def sy_list(N):
    sy_list_flag = []
    for n in range(N):
        op_list = []
        for m in range(N):
            op_list.append(si)
        op_list[n] = sy
        sy_list_flag.append(tensor(op_list))
    return(sy_list_flag)

def sz_list(N):
    sz_list_flag = []
    for n in range(N):
        op_list = []
        for m in range(N):
            op_list.append(si)
        op_list[n] = sz
        sz_list_flag.append(tensor(op_list))
    return(sz_list_flag) 