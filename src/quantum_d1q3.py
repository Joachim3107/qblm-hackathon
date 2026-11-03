"""QLBM D1Q3 (ordres 0 et 1) en encodage d'amplitude, avec Qiskit.

Notations : f^n_{a,i} = population de direction a au nœud i au temps n.
Le champ classique est un tableau f de forme (Q, N) = (3, 2**n).

Registres (Qiskit est little-endian : qubit 0 = bit de poids faible)
  qubits 0 .. n-1   : position i            (N = 2**n nœuds)
  qubits n, n+1     : vitesse a, codée sur 2 bits
                      |00> = repos (c=0), |01> = droite (c=+1), |10> = gauche (c=-1)
                      |11> = inutilisé (amplitude toujours nulle)
  qubit  n+2        : ancilla de la collision (post-sélection sur |0>)
Indice dans le vecteur d'état (sans ancilla) : k = i + N * v
"""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.quantum_info import Statevector

C_VEL = np.array([0, 1, -1])            # c_a
W = np.array([2/3, 1/6, 1/6])           # w_a
CS2 = 1/3
CODE = [0b00, 0b01, 0b10]               # direction a -> état du registre vitesse


# ---------------------------------------------------------------- classique
def collision_matrix(tau, u=0.0):
    """C_ab = (1-w) delta_ab + w * M_ab,  M_ab = w_a (1 + c_a u / cs^2)  (ordre 0 si u=0)."""
    omega = 1 / tau
    M = np.outer(W * (1 + C_VEL * u / CS2), np.ones(3))
    return (1 - omega) * np.eye(3) + omega * M


def classical_step(f, tau, u=0.0):
    fstar = collision_matrix(tau, u) @ f
    return np.array([np.roll(fstar[a], C_VEL[a]) for a in range(3)])


# ---------------------------------------------------------------- encodage
def encode(f):
    """Tableau (3, N) -> vecteur d'état normalisé de 2**(n+2) amplitudes."""
    N = f.shape[1]
    psi = np.zeros(4 * N)
    for a in range(3):
        psi[CODE[a] * N:(CODE[a] + 1) * N] = f[a]
    return psi / np.linalg.norm(psi)


def decode(psi, N, mass):
    """Vecteur d'état -> tableau (3, N), remis à l'échelle par la masse totale (conservée)."""
    f = np.array([np.real(psi[CODE[a] * N:(CODE[a] + 1) * N]) for a in range(3)])
    return f * mass / f.sum()


# ---------------------------------------------------------------- collision
def collision_block_encoding(tau, u=0.0):
    """Porte unitaire 3 qubits (2 vitesse + 1 ancilla) dont le bloc |0>_anc vaut C/s.

    C (3x3) est complétée en C4 (4x4) sur le registre vitesse (état |11> -> 0),
    puis plongée dans l'unitaire U = [[A, sqrt(I - A A^T)], [sqrt(I - A^T A), -A^T]],
    avec A = C4 / s et s = norme de C4 (il faut ||A|| <= 1 pour que U existe).
    """
    C4 = np.zeros((4, 4))
    idx = CODE
    for a in range(3):
        for b in range(3):
            C4[idx[a], idx[b]] = collision_matrix(tau, u)[a, b]
    s = np.linalg.norm(C4, 2) * (1 + 1e-12)
    A = C4 / s

    def sqrtm_psd(X):
        vals, vecs = np.linalg.eigh((X + X.T) / 2)
        return vecs @ np.diag(np.sqrt(np.clip(vals, 0, None))) @ vecs.T

    I4 = np.eye(4)
    U = np.block([[A, sqrtm_psd(I4 - A @ A.T)],
                  [sqrtm_psd(I4 - A.T @ A), -A.T]])
    assert np.allclose(U @ U.T, np.eye(8), atol=1e-10), "U n'est pas unitaire"
    # ordre des qubits de la porte : [vit0, vit1, ancilla] ; ancilla = bit de poids fort
    return UnitaryGate(U, label="collision"), s


# ---------------------------------------------------------------- streaming
def increment(n):
    """|i> -> |i+1 mod 2**n>."""
    qc = QuantumCircuit(n, name="+1")
    for k in reversed(range(1, n)):
        qc.mcx(list(range(k)), k)
    qc.x(0)
    return qc.to_gate()


def streaming_circuit(n):
    """Décale le registre position de c_a, contrôlé par le registre vitesse."""
    inc = increment(n)
    pos, vel = list(range(n)), [n, n + 1]
    qc = QuantumCircuit(n + 2, name="streaming")
    qc.append(inc.control(2, ctrl_state=0b01), vel + pos)            # a = +  -> i+1
    qc.append(inc.inverse().control(2, ctrl_state=0b10), vel + pos)  # a = -  -> i-1
    return qc


# ---------------------------------------------------------------- un pas
def step_circuit(n, tau, u=0.0):
    coll, s = collision_block_encoding(tau, u)
    qc = QuantumCircuit(n + 3)
    qc.append(coll, [n, n + 1, n + 2])
    qc.compose(streaming_circuit(n), qubits=list(range(n + 2)), inplace=True)
    return qc, s


def quantum_step(psi_sys, n, tau, u=0.0):
    """Applique un pas, post-sélectionne l'ancilla sur |0>.
    Renvoie (nouvel état normalisé, probabilité de succès)."""
    qc, _ = step_circuit(n, tau, u)
    dim = 2 ** (n + 2)
    full = np.zeros(2 * dim, dtype=complex)
    full[:dim] = psi_sys                         # ancilla dans |0>
    out = Statevector(full).evolve(qc).data
    kept = out[:dim]                             # composante ancilla = |0>
    p = float(np.vdot(kept, kept).real)
    return kept / np.sqrt(p), p


def run(f0, tau, T, u=0.0):
    N = f0.shape[1]
    n = int(np.log2(N))
    mass = f0.sum()
    psi = encode(f0)
    hist, probs = [f0.copy()], []
    for _ in range(T):
        psi, p = quantum_step(psi, n, tau, u)
        probs.append(p)
        hist.append(decode(psi, N, mass))
    return hist, probs
