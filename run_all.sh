#!/usr/bin/env bash
# Reproduce the maintained deterministic claims in mathematical order.
set -euo pipefail
export PYTHONHASHSEED=0

# Select an interpreter that carries the pinned dependencies, not merely one
# that starts. An explicit choice may be supplied as PY=/path/to/python.
# On systems where `bash` resolves to the Windows Subsystem for Linux, the
# default interpreter is a separate installation and may lack them.
PY="${PY:-python3}"
"$PY" -c "import flint" 2>/dev/null || PY=python
"$PY" -c "import flint" 2>/dev/null || {
  echo "no interpreter with the pinned dependencies found (python-flint missing)."
  echo "install requirements-lock.txt, or run: PY=/path/to/python bash run_all.sh"
  exit 1
}
"$PY" -c "import sys, flint; print('Interpreter:', sys.executable, '| flint', flint.__version__)"

printf '%s\n' '== Canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' '== Decimal interval negative tests =='
"$PY" code/test_decimal_interval.py

printf '%s\n' '== Focused fail-closed regressions =='
"$PY" code/test_fail_closed.py

printf '%s\n' '== Exact geometry and theorem package =='
"$PY" code/verify_theorem_package.py

printf '%s\n' '== Analytic, winding, and chamber package =='
"$PY" code/verify_structural_package.py

printf '%s\n' '== Separate internal analytic determinant reconstruction =='
"$PY" code/check_analytic_det.py

printf '%s\n' '== Separate internal witness reconstruction =='
"$PY" code/check_witness_core.py

printf '%s\n' '== Separate internal exact polygon reconstruction =='
"$PY" code/polygon_indep.py

printf '%s\n' '== Final canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' 'VERDICT: ALL MAINTAINED CHECKS PASSED'
