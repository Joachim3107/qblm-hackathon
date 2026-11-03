"""QLBM D1Q2 en encodage d'amplitude (Qiskit).

Registres (little-endian Qiskit) :
  qubits 0..n-1 : position x (q0 = bit de poids faible)
  qubit n       : vitesse v (0 = f+, 1 = f-)
  qubits n+1, n+2 : ancillas de la collision LCU
Indice du vecteur d'état (sans ancillas) : v * 2**n + x  -> psi = concat(f+, f-).
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate, StatePreparation
from qiskit.quantum_info import Statevector

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])


# ---------- Streaming ----------
def increment(n):
    """|x> -> |x+1 mod 2^n>."""
    qc = QuantumCircuit(n, name="+1")
    for k in reversed(range(1, n)):
        qc.mcx(list(range(k)), k)
    qc.x(0)
    return qc


def streaming_circuit(n):
    inc = increment(n).to_gate()
    qc = QuantumCircuit(n + 1, name="stream")
    pos = list(range(n))
    qc.append(inc.control(1, ctrl_state=0), [n] + pos)            # f+ -> droite
    qc.append(inc.inverse().control(1, ctrl_state=1), [n] + pos)  # f- -> gauche
    return qc


# ---------- Collision (LCU) ----------
def collision_matrix(tau, u):
    w = 1 / tau
    a, b = (1 + u) / 2, (1 - u) / 2
    M = np.array([[a, a], [b, b]])
    return (1 - w) * I2 + w * M


def lcu_terms(tau, u):
    """C = (1-w/2) I + (w/2) X + (wu/2) Z + (wu/2) ZX  ->  [(alpha_k >= 0, U_k unitaire)]."""
    w = 1 / tau
    coeffs = [1 - w / 2, w / 2, w * u / 2, w * u / 2]
    unitaries = [I2, X, Z, Z @ X]
    return [(abs(c), np.sign(c) * U if c != 0 else U) for c, U in zip(coeffs, unitaries)]


def collision_lcu_circuit(n, tau, u):
    """Circuit sur n+3 qubits : agit sur le qubit vitesse (donc sur tous les noeuds à la fois)."""
    terms = lcu_terms(tau, u)
    alphas = np.array([a for a, _ in terms])
    lam = alphas.sum()
    prep = StatePreparation(np.sqrt(alphas / lam))
    anc = [n + 1, n + 2]
    qc = QuantumCircuit(n + 3, name="collide")
    qc.append(prep, anc)
    for k, (_, U) in enumerate(terms):
        qc.append(UnitaryGate(U).control(2, ctrl_state=k), anc + [n])
    qc.append(prep.inverse(), anc)
    return qc, lam


# ---------- Encodage / décodage ----------
def encode(f):
    psi = np.concatenate([f[0], f[1]]).astype(complex)
    return psi / np.linalg.norm(psi)


def decode(psi_sys, mass):
    """Amplitudes -> f, en dénormalisant par la masse totale (conservée)."""
    amp = np.real(psi_sys)
    N = len(amp) // 2
    f = np.array([amp[:N], amp[N:]])
    return f * mass / f.sum()


# ---------- Un pas avec post-sélection (simulation statevector) ----------
def quantum_step(psi_sys, n, tau, u):
    """Retourne (nouvel état normalisé, probabilité de succès de la post-sélection)."""
    coll, _ = collision_lcu_circuit(n, tau, u)
    full = QuantumCircuit(n + 3)
    full.compose(coll, inplace=True)
    full.compose(streaming_circuit(n), qubits=list(range(n + 1)), inplace=True)
    dim = 2 ** (n + 1)
    psi_full = np.zeros(4 * dim, dtype=complex)
    psi_full[:dim] = psi_sys                       # ancillas dans |00>
    out = Statevector(psi_full).evolve(full).data
    kept = out[:dim]                               # post-sélection ancillas = |00>
    p = float(np.vdot(kept, kept).real)
    return kept / np.sqrt(p), p


def run_quantum(f0, tau, u, T):
    n = int(np.log2(f0.shape[1]))
    mass = f0.sum()
    psi = encode(f0)
    history, probs = [f0.copy()], []
    for _ in range(T):
        psi, p = quantum_step(psi, n, tau, u)
        probs.append(p)
        history.append(decode(psi, mass))
    return history, probs
