#!/usr/bin/env bash
# Installe l'environnement du projet dans .venv (n'affecte pas le reste de ton Mac).
set -e
cd "$(dirname "$0")"

PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3,10) else 1)'; then
    PY="$c"; break
  fi
done
if [ -z "$PY" ]; then
  echo "Python >= 3.10 introuvable. Installe-le avec :  brew install python@3.12  puis relance ce script."
  exit 1
fi
echo "Python utilisé : $($PY --version)"

$PY -m venv .venv
./.venv/bin/pip install --upgrade pip -q
./.venv/bin/pip install -r requirements.txt
echo
./.venv/bin/python check_install.py
echo
echo "Pour travailler :  source .venv/bin/activate"
