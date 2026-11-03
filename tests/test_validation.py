import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from qiskit.quantum_info import Statevector
from src import classical_lbm as cl
from src import quantum_lbm as ql


def test_increment():
    n = 3
    for x in range(8):
        sv = Statevector.from_int(x, 2 ** n).evolve(ql.increment(n))
        assert np.argmax(np.abs(sv.data)) == (x + 1) % 8


def test_lcu_decomposition():
    tau, u = 0.8, 0.3
    C = sum(a * U for a, U in ql.lcu_terms(tau, u))
    assert np.allclose(C, ql.collision_matrix(tau, u))


def test_quantum_matches_classical():
    N, tau, u, T = 8, 0.8, 0.3, 10
    f0 = cl.gaussian_init(N, center=2, width=1.5, u=u)
    hc = cl.run(f0, tau, u, T)
    hq, probs = ql.run_quantum(f0, tau, u, T)
    for fc, fq in zip(hc, hq):
        assert np.allclose(fc, fq, atol=1e-10)
    assert all(0 < p <= 1 for p in probs)
