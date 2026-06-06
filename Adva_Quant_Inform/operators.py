

import numpy as np  
import numpy.linalg as LA
import scipy.linalg as sciLA
from qutip import basis 
import scipy.sparse as sparse
from ipywidgets import interactive, interact
from ipywidgets import FloatSlider 
from qutip import tensor, qeye, sigmax, sigmay, sigmaz 
from numpy import (array, pi, cos, sin, ones, size, sqrt, real, mod, append, arange, exp)

X = np.array([[0, 1],
              [1, 0]])

Y = np.array([[0, -1j],
              [1j, 0]])

Z = np.array([[1, 0],
              [0, -1]])

H = np.array([[1, 1],
              [1, -1]]) / np.sqrt(2)

S = np.array([[1, 0],
              [0, 1j]])

T = np.array([[1, 0],
              [0, np.exp(1j*np.pi/4)]])

P0 = np.array([[1, 0],
               [0, 0]])

P1 = np.array([[0, 0],
               [0, 1]])

I = np.identity(2) 

ket0 = 1/np.sqrt(2)*basis(2,0) + 1/np.sqrt(2)*basis(2,1) 
def psi0(N):
    psi0_flag = tensor([ket0 for n in range(N)])
    return(psi0_flag)
    
# helper function for generating basis vectors
def basisvec(n, k):
    v = np.zeros(2**n)
    v[k] = 1
    return v 
    
def rotation(ax,theta):
    return sciLA.expm(-1j * theta/2 * (ax[0]*X + ax[1]*Y + ax[2]*Z)) 
    

def buildSparseGateSingle(n, i, gate):
    """
    Construct a single-qubit gate acting on qubit i in an n-qubit system
    using sparse Kronecker products.

    This embeds a 2×2 quantum gate into the full 2^n-dimensional Hilbert space:

        I ⊗ I ⊗ G ⊗ I ⊗ ... ⊗ I

    where the gate G is applied at position i.

    Parameters
    ----------
    n : Total number of qubits in the system.
    i : Target qubit index (0 ≤ i < n).
    gate :  2×2 quantum gate to apply (e.g., X, H, Z).
    
    Notes
    -----
    - Uses Kronecker products to embed single-qubit operations.
    """

    sgate = sparse.csr_matrix(gate)
    return sparse.kron(
        sparse.kron(
            sparse.identity(2**i, format="csr"),
            sgate
        ),
        sparse.identity(2**(n - i - 1), format="csr")
    )


def buildSparseCNOT(n, ic, it):
    """
    Construct an n-qubit controlled-NOT (CNOT) gate using sparse matrices.
    The CNOT gate is built using projector decomposition:
        CNOT = |0><0|_c ⊗ I_t + |1><1|_c ⊗ X_t

    where:
        - n  : total number of qubits
        - ic : control qubit index (0 ≤ ic < n)
        - it : target qubit index (0 ≤ it < n, it ≠ ic)
    Notes
    -----
    - Uses tensor-product construction with sparse matrices.
    - Computational cost grows exponentially with n (O(2^n)).
    """

    P0_term = buildSparseGateSingle(n, ic, P0)
    P1_term = buildSparseGateSingle(n, ic, P1)
    X_term  = buildSparseGateSingle(n, it, X)
    return P0_term + P1_term @ X_term
    
def U_N_qubits(ops):
    """
    Constructs an N-qubit operator using tensor products.
    """
    U = ops[0]
    for op in ops[1:]:
        U = np.kron(U, op)
    return U

def U_one_gate(V, i, N):
    """
    Applies a single-qubit gate to qubit i in an N-qubit system.

    Parameters
   ...........
    V : Single-qubit gate.
    i : Target qubit index.
    N : Total number of qubits.
    """
    ops = [I] * N
    ops[i] = V
    return U_N_qubits(ops)
    
def controlled_gate(U, control, target, N):
    """
    Controlled-U gate on an N-qubit register.
    Parameters
    ...........
   U:Single-qubit gate
   N: total number of qubits 
    """
    if control == target:
        raise ValueError("Control and target must be different.")

    # Operator acting on the subspace where control qubit is |0⟩
    P0_ops = [
        P0 if i == control else I
        for i in range(N)
    ]

    # Operator acting on the subspace where control qubit is |1⟩
    P1_ops = [
        P1 if i == control else U if i == target else I
        for i in range(N)
    ]

    return U_N_qubits(P0_ops) + U_N_qubits(P1_ops)

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