"""QLBM D1Q3 (ordres 0 et 1) en encodage d'amplitude, avec MyQLM.

Notations : f^n_{a,i} = population de direction a au nœud i au temps n,
champ classique = tableau f de forme (3, N), N = 2**n.

ATTENTION : MyQLM est big-endian (le premier qubit alloué est le bit de POIDS FORT),
contrairement à Qiskit.

Registres, dans l'ordre d'allocation :
  anc  : 1 qubit  (ancilla de collision, post-sélection sur |0>)   <- poids fort
  vel  : 2 qubits (direction a : |00>=repos, |01>=droite, |10>=gauche, |11> inutilisé)
  pos  : n qubits (position i, pos[0] = bit de poids fort)
Indice d'une amplitude : k = anc * 2**(n+2) + v * N + i
Le bloc anc = 0 correspond donc aux 4N premières amplitudes, rangées comme le tableau f.
"""
import numpy as np
from qat.lang.AQASM import Program, QRoutine, AbstractGate, X
from qat.qpus import PyLinalg

C_VEL = np.array([0, 1, -1])            # c_a
W = np.array([2/3, 1/6, 1/6])           # w_a
CS2 = 1/3
CODE = [0b00, 0b01, 0b10]               # direction a -> valeur v du registre vitesse
QPU = PyLinalg()


# ---------------------------------------------------------------- classique (référence)
def collision_matrix(tau, u=0.0):
    """C = (1 - 1/tau) I + M/tau,  M_ab = w_a (1 + c_a u / cs^2)."""
    omega = 1 / tau
    M = np.outer(W * (1 + C_VEL * u / CS2), np.ones(3))
    return (1 - omega) * np.eye(3) + omega * M


def classical_step(f, tau, u=0.0):
    fstar = collision_matrix(tau, u) @ f
    return np.array([np.roll(fstar[a], C_VEL[a]) for a in range(3)])


# ---------------------------------------------------------------- encodage / décodage
def encode(f):
    N = f.shape[1]
    psi = np.zeros(4 * N)
    for a in range(3):
        psi[CODE[a] * N:(CODE[a] + 1) * N] = f[a]
    return psi / np.linalg.norm(psi)


def decode(psi, N, mass):
    f = np.array([np.real(psi[CODE[a] * N:(CODE[a] + 1) * N]) for a in range(3)])
    return f * mass / f.sum()


def state_prep_matrix(psi):
    """Unitaire réelle dont la 1re colonne est psi (réflexion de Householder) : |0> -> |psi>."""
    d = len(psi)
    e0 = np.zeros(d); e0[0] = 1.0
    v = e0 - psi
    if np.linalg.norm(v) < 1e-14:
        return np.eye(d)
    v = v / np.linalg.norm(v)
    return np.eye(d) - 2 * np.outer(v, v)


# ---------------------------------------------------------------- collision (block-encoding)
def collision_unitary(tau, u=0.0):
    """U (8x8) sur [anc, vel0, vel1] dont le bloc anc=0 vaut C4/s (dilatation unitaire)."""
    C4 = np.zeros((4, 4))
    C = collision_matrix(tau, u)
    for a in range(3):
        for b in range(3):
            C4[CODE[a], CODE[b]] = C[a, b]
    s = np.linalg.norm(C4, 2) * (1 + 1e-12)
    A = C4 / s

    def sqrtm_sym(Y):
        vals, vecs = np.linalg.eigh((Y + Y.T) / 2)
        return vecs @ np.diag(np.sqrt(np.clip(vals, 0, None))) @ vecs.T

    I4 = np.eye(4)
    U = np.block([[A, sqrtm_sym(I4 - A @ A.T)],
                  [sqrtm_sym(I4 - A.T @ A), -A.T]])
    assert np.allclose(U @ U.T, np.eye(8), atol=1e-10)
    return U, s


# ---------------------------------------------------------------- streaming
def increment(n):
    """QRoutine |i> -> |i+1 mod 2**n> sur n qubits (big-endian : wire 0 = poids fort)."""
    r = QRoutine()
    w = r.new_wires(n)
    for k in range(n - 1):                 # du bit de poids fort vers le faible
        controls = [w[j] for j in range(k + 1, n)]
        r.apply(X.ctrl(len(controls)), *controls, w[k])
    r.apply(X, w[n - 1])
    return r


# ---------------------------------------------------------------- un pas
def step_program(psi_sys, n, tau, u=0.0):
    """Programme MyQLM : préparation de psi, collision, streaming."""
    prog = Program()
    anc = prog.qalloc(1)
    vel = prog.qalloc(2)
    pos = prog.qalloc(n)

    PREP = AbstractGate("PREP", [], arity=n + 2,
                        matrix_generator=lambda: state_prep_matrix(psi_sys))
    prog.apply(PREP(), vel, pos)

    U, _ = collision_unitary(tau, u)
    COLL = AbstractGate("COLLISION", [], arity=3, matrix_generator=lambda: U)
    prog.apply(COLL(), anc, vel)

    inc = increment(n)
    # direction +  (vel = |01>) : i -> i+1      (contrôle sur vel0 = 0 : on l'entoure de X)
    prog.apply(X, vel[0])
    prog.apply(inc.ctrl(2), vel[0], vel[1], pos)
    prog.apply(X, vel[0])
    # direction -  (vel = |10>) : i -> i-1
    prog.apply(X, vel[1])
    prog.apply(inc.dag().ctrl(2), vel[0], vel[1], pos)
    prog.apply(X, vel[1])
    return prog


def statevector(prog):
    circ = prog.to_circ()
    res = QPU.submit(circ.to_job())
    vec = np.zeros(2 ** circ.nbqbits, dtype=complex)
    for sample in res:
        vec[sample.state.int] = sample.amplitude
    return vec


def quantum_step(psi_sys, n, tau, u=0.0):
    """Un pas + post-sélection anc = 0. Renvoie (nouvel état normalisé, proba de succès)."""
    vec = statevector(step_program(psi_sys, n, tau, u))
    kept = vec[:2 ** (n + 2)]              # bloc anc = 0
    p = float(np.vdot(kept, kept).real)
    return np.real(kept) / np.sqrt(p), p


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
