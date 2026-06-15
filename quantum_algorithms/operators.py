

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
I8 = np.eye(8, dtype=complex)

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


def projectors(dim):
    """
    Generate computational basis projectors {|i><i|} with the given dimension.
    """
    projectors = []
    for i in range(dim):
        ket = np.zeros(dim, dtype=complex)
        ket[i] = 1
        P = np.outer(ket, ket)
        projectors.append(P)
    return projectors
    
def rotation_gate(theta, n):
    """
    This function implements a unitary rotation of a single qubit
    by an angle `theta` around an axis `n` on the Bloch sphere.

    The rotation generator is constructed as N = n · σ,
    where σ = (X, Y, Z) are the Pauli matrices.

    Parameters
    theta : Rotation angle
    n : Rotation axis
    """
    nx, ny, nz = n
    N = nx * X + ny * Y + nz * Z
    R = np.cos(theta / 2) * I - 1j * np.sin(theta / 2) * N
    return R
    
def U_N_qubits(ops):
    """
    Constructs an N-qubit operator using tensor products.

    Parameters
    ops : single-qubit operators.
    """
    U = ops[0]
    for op in ops[1:]:
        U = np.kron(U, op)
    return U


def U_one_gate(V, i, N):
    """
    Applies a single-qubit gate to qubit i
    in an N-qubit system.

    Parameters
    V : Single-qubit gate.
    i : Target qubit index.
    N : Total number of qubits.
    """
    ops = [I] * N
    ops[i] = V
    return U_N_qubits(ops)


def U_two_gates(V, W, i, j, N):
    """
    Applies two single-qubit gates to an N-qubit system.
    If i != j:
        applies V on qubit i and W on qubit j.
    If i == j:
        applies the composed gate V @ W on qubit i,
        preserving operator ordering.
    """
    ops = [I] * N
    if i == j:
        ops[i] = V @ W
    else:
        ops[i] = V
        ops[j] = W
    return U_N_qubits(ops)

def rho(states, probabilities):
    """
    Constructs a density matrix from pure states.
    Parameters
    states  and probabilities 
    """
    return sum(p * np.outer(psi, psi.conj())
               for psi, p in zip(states, probabilities))

def evolve(state, U):
    """
    Evolves a quantum state using a unitary operator.
    Parameters
    state : State vector or density matrix.
    U :  Unitary operator.
    """
    if state.ndim == 1:
        # Pure state evolution
        return U @ state
    elif state.ndim == 2:
        # Density matrix evolution
        return U @ state @ U.conj().T
    else:
        raise ValueError("State must be a vector or a density matrix")


def controlled_gate(U, control, target, N):
    """
    Controlled-U gate on an N-qubit register.
    Implements the projector decomposition:
        C_U = P0(control) ⊗ I  +  P1(control) ⊗ U(target)
    """
    if control == target:
        raise ValueError("Control and target must be different")
    # Operator acting on the subspace where the control qubit is |0⟩
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

def normalize_state(psi):
    """
    Normalize a pure state vector |psi>.
    """
    norm = np.linalg.norm(psi)
    if np.isclose(norm, 0):
        raise ValueError("State vector has zero norm.")
    return psi/norm 

def born_rule_probs(rho, projectors):
    """
    Compute measurement outcome probabilities using the Born rule:
    p(i)= Tr(Pi * rho) for each projector
    """
    probs = np.array([np.real(np.trace(Pi @ rho)) for Pi in projectors])
    probs = np.clip(probs, 0, 1)
    probs = probs / np.sum(probs)
    return probs


def sample_from_probs(probs):
    """
    Return a sampled index based on probs.
    """
    return np.random.choice(len(probs), p=probs)

def measure_pure_state(psi, projectors):
    """
    Measure pure state |psi> using projectors.
    Returns:
        outcome (int)
        psi_post (np.ndarray)
        probs (np.ndarray)
    """
    psi = normalize_state(psi)
    probs = born_probs_pure(psi, projectors)
    outcome = sample_from_probs(probs)
    Pk = projectors[outcome]
    psi_post_unnormalized = Pk @ psi
    norm_post = np.linalg.norm(psi_post_unnormalized)
    if np.isclose(norm_post, 0):
        raise ValueError("Outcome probability ~0 (numerical issue).")
    psi_post = psi_post_unnormalized / norm_post
    
    return outcome, psi_post, probs
    
def measurement_density_matrix(rho, projectors):
    """
    Perform measurement using the given projectors
    """
    probs = born_rule_probs(rho, projectors)
    outcome = sample_from_probs(probs)
    Pk = projectors[outcome]
    numerator = Pk @ rho @ Pk
    denom = np.trace(numerator)
    if np.isclose(denom, 0):
        raise ValueError("Outcome probability ~0 (numerical issue).")
    rho_post = numerator / denom
    return outcome, rho_post, probs

def initial_state(n):
    """
    Prepared the initial state
    """
    total = n + 1
    state = np.zeros(2**total, dtype=complex)
    state[1] = 1.0  # basis index where ancilla=1 and inputs all zero
    return state

def apply_hadamards(state, total_qubits):
    """
    Apply Hadamards gate on the prepared initial state to create a superposition.
    """
    H_full = U_N_qubits([H] * total_qubits)
    return H_full.dot(state)

def sample_probs(probs, shots, rng=None):
    """
    Sample measurement outcomes based on a given probability distribution.
    
    It draws a specified number of random samples ("shots") according to the
    probability distribution `probs`,

    """
    if rng is None:
        rng = np.random.default_rng()
    outcomes = rng.choice(len(probs), size=shots, p=probs)
    return Counter(outcomes)
    
def oracle_function(f, n):
    """
    Build a function that applies the oracle operator U_f to the statevector of n+1 qubits.
    The oracle implements the transformation:
        U_f |x⟩|y⟩ = |x⟩|y ⊕ f(x)⟩
    Parameters
        f : Boolean function f(x) -> {0, 1} for x in [0, 2^n)
        n : Number of input qubits
    """

    def apply_Uf(state):
        new = np.copy(state)
        for x in range(2**n):
            fx = f(x)
            idx0 = (x << 1) | 0
            idx1 = (x << 1) | 1
            if fx == 1:
                # swap amplitudes between ancilla 0 and 1 for this x
                new[idx0], new[idx1] = state[idx1], state[idx0]
        return new   
    return apply_Uf  

def f_constant_0(x):
    return 0 

def f_constant_1(x):
    return 1

def f_balanced_parity(x):
    return x % 2  # 0 for even, 1 for odd 

def measure_probs_first_n(state, n):
    """Compute prob distribution over first n qubits (sum over ancilla)."""
    probs = np.zeros(2**n)
    for x in range(2**n):
        # Apply bitwise operations to find the correct index for each state
        idx0 = (x << 1) | 0  # ancilla = 0
        idx1 = (x << 1) | 1  # ancilla = 1
        probs[x] = np.abs(state[idx0])**2 + np.abs(state[idx1])**2
    return probs 

def sample_measurements_input(state, n, shots, rng=None):
    """
    Measurement outcomes from the full-register distribution given by state,
    then aggregate counts over the input register (i.e., ignore ancilla).
    """
    if rng is None:
        rng = np.random.default_rng()
    probs_full = np.abs(state)**2
    probs_full = probs_full / probs_full.sum()
    samples = rng.choice(len(probs_full), size=shots, p=probs_full)
    input_samples = samples >> 1   # removes ancilla qubit (shift right)
    return Counter(input_samples) 

def bloch_vector(rho):
    """
    Compute the Bloch vector (rX, rY, rZ) for a single-qubit density matrix rho.
    r_J = Tr(rho * J), J = X, Y, Z
    """
    rX = np.real(np.trace(rho @ X))
    rY = np.real(np.trace(rho @ Y))
    rZ = np.real(np.trace(rho @ Z))
    return np.array([rX, rY, rZ])
    

def rotation_channel(p, R):
    """
    Random unitary single-qubit channel using rotation R.
    Returns list of Kraus operators [M0, M1].
    """
    M0 = np.sqrt(1-p) * I
    M1 = np.sqrt(p) * R
    return [M0, M1]


def ket0():
    return np.array([1, 0], dtype=complex)

def ket1():
    return np.array([0, 1], dtype=complex)

def ket_plus():
    return (ket0() + ket1()) / np.sqrt(2)

def ket_minus():
    return (ket0() - ket1()) / np.sqrt(2)

def dm(psi):
    """
    Construct a density matrix from a pure state |psi⟩.

    ρ = |psi⟩⟨psi|
    """
    psi = psi/np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def random_pure_state(rng=None):
    """
    Generate a random single-qubit pure state |psi⟩.
    Useful for Monte-Carlo simulations.
    """
    rng = np.random.default_rng() if rng is None else rng
    v = rng.normal(size=2) + 1j * rng.normal(size=2)
    v = v / np.linalg.norm(v)
    return v


def E1_rho(psi, p):
    """
    Q3: One-qubit bit flip channel:
        E1(rho) = (1-p)rho + pXrhoX
    """
    rho = dm(psi)
    return (1 - p) * rho + p * (X @ rho @ X)

def deutsch_jozsa(n, f):
    """
    Deutsch–Jozsa Algorithm(DJA) the Boolean function 
        is constant or balanced using a single oracle query. 
    Parameters:
        n : Number of input qubits
        f : Oracle function
    """
    total_qubits = n + 1
    state = initial_state(n)
    state = apply_hadamards(state, total_qubits)
    U = oracle_function(f, n)
    state = U(state)
    H_first_n = U_N_qubits([H] * n + [I])
    state = H_first_n @ state
    return state  

def deutsch_jozsa_error1(n, f, theta, target_qubit, axis):
    """
    DJA with a single-qubit rotation error applied before the first Hadamard gates.
    """
    total_qubits = n + 1
    state = initial_state(n)
    R = rotation_gate(theta, axis)
    state = U_one_gate(R, target_qubit, total_qubits) @ state
    state = apply_hadamards(state, total_qubits)
    U = oracle_function(f, n)
    state = U(state)
    H_first_n = U_N_qubits([H]*n + [I])
    state = H_first_n @ state

    return state

def deutsch_jozsa_error2(n, f, theta, target_qubit, axis):
    """
    DJA with a single-qubit rotation error applied after the first Hadamard gates.
    """
    total_qubits = n + 1
    state = initial_state(n)
    state = apply_hadamards(state, total_qubits)
    R = rotation_gate(theta, axis)
    state = U_one_gate(R, target_qubit, total_qubits) @ state
    U = oracle_function(f, n)
    state = U(state)
    H_first_n = U_N_qubits([H]*n + [I])
    state = H_first_n @ state
    
    return state

def deutsch_jozsa_error3(n, f, theta, target_qubit, axis):
    """
    DJA with a single-qubit error applied after the oracle U_f.
    """
    total_qubits = n + 1
    state = initial_state(n)
    state = apply_hadamards(state, total_qubits)
    U = oracle_function(f, n)
    state = U(state)
    R = rotation_gate(theta, axis)
    state = U_one_gate(R, target_qubit, total_qubits) @ state
    H_first_n = U_N_qubits([H]*n + [I])
    state = H_first_n @ state

    return state

def deutsch_jozsa_error4(n, f, theta, target_qubit, axis):
    """
    Deutsch–Jozsa algorithm with a single-qubit rotation error
    applied after the final Hadamard gates.
    """
    total_qubits = n + 1
    state = initial_state(n)
    state = apply_hadamards(state, total_qubits)
    U = oracle_function(f, n)
    state = U(state)
    H_first_n = U_N_qubits([H]*n + [I])
    state = H_first_n @ state
    R = rotation_gate(theta, axis)
    state = U_one_gate(R, target_qubit, total_qubits) @ state

    return state
    


# Sparse helper functions
def buildSparseGateSingle(n, i, gate):
    return sparse.kron(sparse.kron(sparse.identity(2**i), gate), sparse.identity(2**(n-i-1)))

def buildSparseCNOT(n, ic, it):
    return buildSparseGateSingle(n, ic, P0) + buildSparseGateSingle(n, ic, P1) @ buildSparseGateSingle(n, it, X)

def dm_sparse(psi):
    """Density matrix from state vector (sparse)"""
    return psi @ psi.getH()

def ket0_sparse(n=1):
    """n-qubit |0>"""
    return sparse.csr_matrix(np.array([[1], [0]], dtype=complex)) if n==1 else sparse.kron(ket0_sparse(), sparse.identity(2**(n-1), dtype=complex))

# Multi-qubit bit-flip Kraus (sparse)
def bit_flip_kraus_nqubits_sparse(p, n):
    single_ops = [np.sqrt(1-p) * I, np.sqrt(p) * X]
    # generate all combinations
    kraus_ops = []
    for combo in product(single_ops, repeat=n):
        K = combo[0]
        for E in combo[1:]:
            K = sparse.kron(K, E)
        kraus_ops.append(K)
    return kraus_ops

def buildSparseCNOT(n, ic, it):
    P0ic = buildSparseGateSingle(n, ic, P0)
    P1ic = buildSparseGateSingle(n, ic, P1)
    Xit  = buildSparseGateSingle(n, it, X)
    return P0ic + P1ic @ Xit


# helper function for initializing all qubits in state zero
def initRegisterPsi(n):
    return basisvec(n,0)

def initRegisterRho(n):
    ini = basisvec(n,0)
    return np.outer(ini.conj(),ini)

