# qblm-hackathon

**Quantum Lattice Boltzmann Methods (QLBM) pour la simulation de fluides** — hackathon de 72 h, Team 2.

On simule l'équation de Boltzmann sur réseau (schéma D1Q3, `N = 2^k` nœuds périodiques) **entièrement sur un circuit quantique** (Qiskit) : encodage d'amplitude, collision par *block-encoding*, streaming par incréments contrôlés, post-sélection et reconstruction.

## Résultats principaux

Valeurs issues des notebooks (`k = 3`, donc N = 8 nœuds ; T = 3 pas ; τ = 0,8).

| | Diffusion (FABLE) | Advection-diffusion (FABLE) | Changement de variables + LCU |
|---|---|---|---|
| Qubits | 14 | 14 | **8** |
| Collision | 22 CNOT, 3 ancillas / pas | 22 CNOT, 3 ancillas / pas | **8 CNOT, 1 ancilla / pas** |
| Probabilité de succès (3 pas) | 2,2·10⁻⁴ | 2,4·10⁻⁴ | **0,904** (diff.) / **0,978** (adv.) |
| Circuit exact vs LBM classique | ~10⁻¹⁵ | ~10⁻¹⁴ | ~10⁻¹⁵ |

- **Advection-diffusion** : avec l'équilibre d'ordre 1, `f_eq = ω ρ (1 + c u / c_s²)` et `u = j/ρ`, le produit `ρu = j` rend l'équilibre **linéaire** (`f_eq = E f`). Le circuit est donc le même qu'en diffusion : seuls les angles de FABLE changent.
- **Optimisation** : on pose `G = W^{-1/2} F`. La collision devient alors symétrique et de norme 1 (contre 1,345 avant), et s'écrit `a·I + b·R`, où `R` est une réflexion. Cette combinaison linéaire de deux unitaires (LCU) donne un block-encoding exact, de facteur 1.

## Arborescence

```
├── notebooks/
│   ├── 00_sujet_encadrants/      sujet et notebooks fournis (LBM, fluides, FABLE, Carleman, D1Q3)
│   ├── 01_diffusion_FABLE/       diffusion D1Q3 : circuit complet avec FABLE (référence)
│   ├── 02_advection_diffusion/   advection-diffusion, modèle u = j/ρ
│   ├── 03_comparatif/            diffusion vs advection à N égal (coût, succès, précision)
│   ├── 04_optimisation_LCU/      changement de variables + LCU (méthode principale)
│   ├── outils_myqlm/             lecture pas à pas de src/myqlm_d1q3.py
│   └── archives/                 anciennes versions
├── src/        code réutilisable (LBM classique, circuits Qiskit, version MyQLM)
├── tests/      validation quantique = classique (pytest)
├── figures/    figures exportées des notebooks + schémas de circuits
├── docs/       notes explicatives (PDF)
└── slides/     présentation finale (Manim, HTML)
```

## Ordre de lecture conseillé

1. `notebooks/00_sujet_encadrants/lbm.ipynb` puis `LBM_D1Q3.ipynb` : la LBM classique.
2. `notebooks/00_sujet_encadrants/block_encoding.ipynb` : FABLE.
3. `notebooks/01_diffusion_FABLE/LBM_D1Q3_Circuit_tot.ipynb` : premier circuit complet.
4. `notebooks/02_advection_diffusion/LBM_D1Q3_Advection.ipynb`.
5. `notebooks/03_comparatif/Comparatif_diffusion_advection.ipynb`.
6. `notebooks/04_optimisation_LCU/LBM_D1Q3_LCU.ipynb`.

## Présentation

`slides/LBM_presentation_v10.html` : la présentation finale (animations Manim). Ouvrir le fichier dans un navigateur. Il est volumineux (~14 Mo) : GitHub ne l'affiche pas, il faut le télécharger.

## Installation

```bash
git clone https://github.com/Joachim3107/qblm-hackathon.git
cd qblm-hackathon
./setup.sh                 # crée .venv, installe requirements.txt, lance une démo
source .venv/bin/activate
pytest -q                  # vérifie quantique = classique
jupyter lab
```

Environnement MyQLM (séparé, Python 3.11 à 3.13) :

```bash
./setup_myqlm.sh
source .venv-myqlm/bin/activate
python check_myqlm.py
pytest -q tests/test_myqlm_d1q3.py
```

## Conventions

- Notations D1Q3 : directions `(0, −, +)`, `c = (0, −1, +1)`, `ω = (2/3, 1/6, 1/6)`, `c_s² = 1/3`. La composante `(a, i)` de `F` est à l'indice `a·N + i`.
- Qiskit est *little-endian* : la position est sur `q_0 … q_{k−1}`, la direction sur `q_k, q_{k+1}`, avec `0 → |00⟩`, `− → |01⟩`, `+ → |10⟩` et une direction vide `∅ → |11⟩`.
- MyQLM est *big-endian* : le premier qubit alloué est le bit de poids fort.

## Bibliographie

- FABLE — Camps & Van Beeumen, arXiv:2205.00081
- S-FABLE — Kuklinski & Rempfer, arXiv:2401.04234
- Camps, Lin, Van Beeumen, Yang, arXiv:2203.10236
- Fang, Lin, Tong, *Quantum* 7, 955 (2023)
- Tiwari et al. (Ansys/IonQ), arXiv:2504.10870
- Zeng et al., arXiv:2609.16721
- Budinski, arXiv:2103.03804
- Bediche et al., arXiv:2606.16514
- Holmes & Matsuura, arXiv:2005.04351
