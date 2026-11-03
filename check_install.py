"""Vérifie l'installation et lance une mini-démo QLBM (8 noeuds, 4 qubits + 2 ancillas)."""
import numpy as np, qiskit, qiskit_aer, matplotlib
from src import classical_lbm as cl, quantum_lbm as ql

print(f"qiskit {qiskit.__version__} | qiskit-aer {qiskit_aer.__version__} | numpy {np.__version__}")
N, tau, u, T = 8, 0.8, 0.3, 10
f0 = cl.gaussian_init(N, center=2, width=1.5, u=u)
hc = cl.run(f0, tau, u, T)
hq, probs = ql.run_quantum(f0, tau, u, T)
err = max(np.abs(a - b).max() for a, b in zip(hc, hq))
print(f"Erreur max quantique vs classique sur {T} pas : {err:.2e}")
print(f"Proba de succès par pas : {np.mean(probs):.3f} | cumulée sur {T} pas : {np.prod(probs):.3e}")
print("OK, tout est prêt." if err < 1e-8 else "PROBLÈME : écart anormal.")
