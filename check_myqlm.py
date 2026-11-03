"""Vérifie l'installation MyQLM et lance la démo QLBM D1Q3 (16 nœuds, 7 qubits)."""
import numpy as np
import qat
from qat.lang.AQASM import Program, H, CNOT
from qat.qpus import PyLinalg
from src import myqlm_d1q3 as q

prog = Program(); r = prog.qalloc(2)
prog.apply(H, r[0]); prog.apply(CNOT, r[0], r[1])
bell = {str(s.state): round(s.probability, 3) for s in PyLinalg().submit(prog.to_circ().to_job())}
print("MyQLM OK, état de Bell :", bell)

N, tau, T = 16, 0.8, 8
x = np.arange(N)
f0 = np.outer(q.W, np.ones(1)) * np.exp(-(x - 4) ** 2 / (2 * 1.6 ** 2))
hist, probs = q.run(f0, tau, T)
fc = f0.copy()
for _ in range(T):
    fc = q.classical_step(fc, tau)
err = np.abs(fc - hist[-1]).max()
print(f"QLBM D1Q3 ({T} pas) : erreur max quantique vs classique = {err:.1e}")
print(f"Proba de succès par pas ≈ {np.mean(probs):.3f}, cumulée = {np.prod(probs):.2e}")
print("OK, tout est prêt." if err < 1e-8 else "PROBLÈME : écart anormal.")
