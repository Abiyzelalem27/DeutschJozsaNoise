



import numpy as np
from collections import Counter
import itertools 
from scipy import sparse 
import scipy
import matplotlib.pyplot as plt 
import math
import random 

I = np.array([[1, 0],
              [0, 1]], dtype=complex)
I8 = np.eye(8, dtype=complex)
X = np.array([[0, 1],
              [1, 0]], dtype=complex)
Y = np.array([[0, -1j],
              [1j,  0]], dtype=complex)
Z = np.array([[1,  0],
              [0, -1]], dtype=complex)
H = 1 / np.sqrt(2) * np.array([[1,  1],
                               [1, -1]], dtype=complex)
P0 = np.array([[1, 0],
               [0, 0]], dtype=complex)
P1 = np.array([[0, 0],
               [0, 1]], dtype=complex)

# 3-qubit identity
I8 = np.kron(np.kron(I,I),I)

def clean_probs(probs):
    probs = np.real(probs)
    probs = np.clip(probs, 0, None)
    probs = probs / np.sum(probs)
    return probs

def normalize(v):
    return v / np.linalg.norm(v)

def rho_from_state(psi):
    return np.outer(psi, np.conjugate(psi))


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
    

def doMeasurement(state, projectors):
    """
     Perform a projective measurement on a quantum state (pure or density matrix)
    using a list of projectors
    Returns:
        outcome :Index of the measured outcome
        post_state : Collapsed state after measurement (normalized)
        probs : Probability 
    """
    # Determine if pure state or density matrix
    pure = state.ndim == 1

    if pure:
        # Normalize pure state
        state = state / np.linalg.norm(state)
        # Probabilities via Born rule
        probs = np.array([np.real(np.vdot(state, P @ state)) for P in projectors])
    else:
        # Density matrix case
        probs = np.array([np.real(np.trace(P @ state)) for P in projectors])

    # Normalize probabilities (safety)
    probs = np.clip(probs, 0, 1)
    probs /= probs.sum()
    # Sample outcome
    outcome = np.random.choice(len(projectors), p=probs)
    Pk = projectors[outcome]

    # Collapse state
    if pure:
        post_unnorm = Pk @ state
        norm_post = np.linalg.norm(post_unnorm)
        if np.isclose(norm_post, 0):
            raise ValueError("Outcome probability ~0 (numerical issue).")
        post_state = post_unnorm / norm_post
    else:
        post_state = Pk @ state @ Pk
        denom = np.trace(post_state)
        if np.isclose(denom, 0):
            raise ValueError("Outcome probability ~0 (numerical issue).")
        post_state = post_state / denom
    return outcome, post_state, probs


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


def apply_kraus(rho, kraus_ops):
    """
    Apply a quantum channel to a density matrix using Kraus operators.
    """
    rho_out = np.zeros_like(rho, dtype=complex)
    for E in kraus_ops:
        rho_out += E @ rho @ E.conj().T  # E rho E†
    
    return rho_out 


def depolarizing_kraus(p): 
    """
    Depolarizing channel.

    Kraus operators:
        E0 = sqrt(1 - 3p/4) I
        E1 = sqrt(p/4) X
        E2 = sqrt(p/4) Y
        E3 = sqrt(p/4) Z
    """
    E0 = np.sqrt(1 - 3*p/4) * I
    E1 = np.sqrt(p/4) * X
    E2 = np.sqrt(p/4) * Y
    E3 = np.sqrt(p/4) * Z
    return [E0, E1, E2, E3]


def dm(psi):
    """
    Construct a density matrix from a pure state |psi⟩.

    ρ = |psi⟩⟨psi|
    """
    psi = psi/np.linalg.norm(psi)
    return np.outer(psi, psi.conj())
    
def measure_probs_first_n_rho(rho, n):
    probs = np.zeros(2**n)

    for x in range(2**n):
        idx0 = (x << 1) | 0
        idx1 = (x << 1) | 1

        probs[x] = np.real(
            rho[idx0, idx0] +
            rho[idx1, idx1]
        )

    return probs

    
def single_qubit_channel_n_register(kraus_single, n, target):
    """
    Lift single-qubit Kraus operators to act on qubit 'target' in an n-qubit register.

    Parameters:
        kraus_single : list of 2x2 Kraus operators for the single qubit
        n: total number of qubits in the register
        target: index of qubit to apply the channel (0-based)

    Returns:
        list of 2^n x 2^n Kraus operators acting on the full register
    """
    kraus_n = []

    for K in kraus_single:
        full_op = np.array([[1]], dtype=complex)  
        for qubit in range(n):
            op = K if qubit == target else I
            full_op = np.kron(full_op, op)
        kraus_n.append(full_op)
    return kraus_n 

def oracle_matrix(f, n):
    total_qubits = n + 1
    dim = 2**total_qubits
    U = np.zeros((dim, dim), dtype=complex)

    for basis in range(dim):
        x = basis >> 1
        y = basis & 1
        y_new = y ^ f(x)
        new_basis = (x << 1) | y_new
        U[new_basis, basis] = 1

    return U

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



def deutsch_jozsa_depolarizing1(n, f, p, target_qubit):
    """
    Depolarizing noise at E1:
    Before the first Hadamard gates.
    """
    total_qubits = n + 1

    rho = dm(initial_state(n))  # initial density matrix

    K = depolarizing_kraus(p)    # depolarizing noise E1
    K_full = single_qubit_channel_n_register(K, total_qubits, target_qubit)
    rho = apply_kraus(rho, K_full)

    H_all = U_N_qubits([H] * total_qubits)
    rho = evolve(rho, H_all) # first Hadamards on all qubits

    U_f = oracle_matrix(f, n)  # oracle
    rho = evolve(rho, U_f)   

    H_first_n = U_N_qubits([H] * n + [I])
    rho = evolve(rho, H_first_n)  # final Hadamards on input qubits only

    return rho
    

def deutsch_jozsa_depolarizing2(n, f, p, target_qubit):
    """
    Depolarizing noise at E2:
    After the first Hadamard gates.
    """
    total_qubits = n + 1

    state = initial_state(n)
    rho = dm(state)

    H_all = U_N_qubits([H] * total_qubits)
    rho = evolve(rho, H_all)

    K = depolarizing_kraus(p)  # depolarizing noise
    K_full = single_qubit_channel_n_register(K, total_qubits, target_qubit )

    rho = apply_kraus(rho, K_full)

    U_f = oracle_matrix(f, n)
    rho = evolve(rho, U_f)

    H_first_n = U_N_qubits([H] * n + [I])
    rho = evolve(rho, H_first_n)
    
    return rho
    

def deutsch_jozsa_depolarizing3(n, f, p, target_qubit):
    """
    Depolarizing noise at E3:
    Immediately after the oracle.
    """
    total_qubits = n + 1

    state = initial_state(n)
    rho = dm(state)

    H_all = U_N_qubits([H] * total_qubits)
    rho = evolve(rho, H_all)

    U_f = oracle_matrix(f, n)
    rho = evolve(rho, U_f)
 
    K = depolarizing_kraus(p) # depolarizing noise
    K_full = single_qubit_channel_n_register(
        K, total_qubits, target_qubit
    )

    rho = apply_kraus(rho, K_full)

    H_first_n = U_N_qubits([H] * n + [I])
    rho = evolve(rho, H_first_n)

    return rho

def deutsch_jozsa_depolarizing4(n, f, p, target_qubit):
    """
    Depolarizing noise at E4:
    After the final Hadamard gates.
    """
    total_qubits = n + 1

    state = initial_state(n)
    rho = dm(state)

    H_all = U_N_qubits([H] * total_qubits)
    rho = evolve(rho, H_all)

    U_f = oracle_matrix(f, n)
    rho = evolve(rho, U_f)

    H_first_n = U_N_qubits([H] * n + [I])
    rho = evolve(rho, H_first_n)

    K = depolarizing_kraus(p) # depolarizing noise
    K_full = single_qubit_channel_n_register(
        K, total_qubits, target_qubit
    )

    rho = apply_kraus(rho, K_full)

    return rho 




