# QLBM Hackathon — Team 2

Quantum Lattice Boltzmann Methods pour la simulation de fluides.
Point de départ : advection-diffusion D1Q2, encodage d'amplitude, streaming par incrémenteur, collision par LCU.

## Installation (une fois)

```bash
git clone <URL_DU_REPO>
cd qlbm-hackathon
./setup.sh
source .venv/bin/activate
```

Le script crée un environnement isolé `.venv`, installe les dépendances et lance une démo de vérification.

## Structure

| Fichier | Rôle |
|---|---|
| `src/classical_lbm.py` | Référence classique numpy (vérité terrain) |
| `src/quantum_lbm.py` | Circuits Qiskit : encodage, streaming, collision LCU, post-sélection |
| `tests/test_validation.py` | Vérifie quantique = classique |
| `check_install.py` | Démo rapide + vérification de l'installation |
| `notebooks/` | Explorations, figures |
| `figures/` | Figures pour le pitch |

## Commandes utiles

```bash
python check_install.py     # démo
pytest -q                   # tests
jupyter lab                 # notebooks
```

## Conventions d'équipe

- Une branche par tâche (`git switch -c streaming-qft`), merge via pull request.
- `main` doit toujours marcher : lancer `pytest -q` avant de merger.
- Qiskit est little-endian : qubit 0 = bit de poids faible.

## Environnement MyQLM

Un second environnement, séparé de Qiskit, pour travailler avec MyQLM (Eviden).
MyQLM 1.14 demande Python 3.11, 3.12 ou 3.13.

```bash
./setup_myqlm.sh                    # crée .venv-myqlm et le kernel Jupyter « Python (MyQLM) »
source .venv-myqlm/bin/activate
python check_myqlm.py               # démo QLBM D1Q3 en MyQLM
pytest -q tests/test_myqlm_d1q3.py  # validation quantique = classique
```

| Fichier | Rôle |
|---|---|
| `src/myqlm_d1q3.py` | QLBM D1Q3 en MyQLM : encodage, collision (block-encoding), streaming |
| `tests/test_myqlm_d1q3.py` | Vérifie MyQLM = LBM classique aux ordres 0 et 1 |

Attention : MyQLM est big-endian (premier qubit alloué = bit de poids fort), Qiskit est little-endian.
