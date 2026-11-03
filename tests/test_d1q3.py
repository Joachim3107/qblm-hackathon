import sys, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from src import quantum_d1q3 as q


def gaussian(N, u=0.0):
    x = np.arange(N)
    rho = np.exp(-(x - N / 4) ** 2 / (2 * (N / 10) ** 2))
    return np.outer(q.W * (1 + q.C_VEL * u / q.CS2), np.ones(1)) * rho


def check(u):
    N, tau, T = 16, 0.8, 12
    f = gaussian(N, u)
    hq, probs = q.run(f, tau, T, u)
    fc = f.copy()
    for t in range(T):
        fc = q.classical_step(fc, tau, u)
        assert np.allclose(fc, hq[t + 1], atol=1e-10), f"écart au pas {t+1}"
    return probs


def test_order0_diffusion():
    check(0.0)


def test_order1_advection_diffusion():
    check(0.1)
