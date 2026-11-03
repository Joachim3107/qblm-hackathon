#!/usr/bin/env bash
# Crée l'environnement MyQLM dans .venv-myqlm (séparé de l'environnement Qiskit .venv).
# MyQLM 1.14 demande Python 3.11, 3.12 ou 3.13 sur Mac Apple Silicon.
set -e
cd "$(dirname "$0")"

PY=""
for c in python3.12 python3.13 python3.11; do
  if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
  echo "Aucun Python 3.11 à 3.13 trouvé."
  echo "Installe-le avec :  brew install python@3.12   puis relance ./setup_myqlm.sh"
  exit 1
fi
echo "Python utilisé : $($PY --version)"

$PY -m venv .venv-myqlm
./.venv-myqlm/bin/pip install --upgrade pip -q
./.venv-myqlm/bin/pip install -r requirements-myqlm.txt
./.venv-myqlm/bin/python -m ipykernel install --user --name qlbm-myqlm --display-name "Python (MyQLM)"
echo
./.venv-myqlm/bin/python check_myqlm.py
echo
echo "Pour travailler :  source .venv-myqlm/bin/activate"
echo "Dans VS Code / Jupyter : choisis le kernel « Python (MyQLM) »"
