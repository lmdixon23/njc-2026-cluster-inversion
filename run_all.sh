#!/usr/bin/env bash
# Reproduce the maintained deterministic claims in mathematical order.
set -euo pipefail
export PYTHONHASHSEED=0

# Select an interpreter that carries the pinned dependencies, not merely one
# that starts. An explicit choice may be supplied as PY=/path/to/python.
# On systems where `bash` resolves to the Windows Subsystem for Linux, the
# default interpreter is a separate installation and may lack them.
if [[ ${PY+x} ]]; then
  # An explicit interpreter is a contract: do not silently substitute another.
  "$PY" -c "import flint" 2>/dev/null || {
    echo "FAIL: explicit PY cannot import python-flint; use the pinned environment."
    exit 1
  }
else
  PY=python3
  "$PY" -c "import flint" 2>/dev/null || PY=python
  "$PY" -c "import flint" 2>/dev/null || {
    echo "FAIL: no interpreter with python-flint; install requirements-lock.txt."
    exit 1
  }
fi
"$PY" -c "import sys; sys.exit('FAIL: verification requires assertions; unset PYTHONOPTIMIZE.') if sys.flags.optimize else None"
"$PY" -c "import sys, flint; print('Interpreter:', sys.executable, '| flint', flint.__version__)"

printf '%s\n' '== Verification runtime regressions =='
"$PY" code/test_verification_runtime.py

printf '%s\n' '== Canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' '== Decimal interval negative tests =='
"$PY" code/test_decimal_interval.py

printf '%s\n' '== Focused fail-closed regressions =='
"$PY" code/test_fail_closed.py
"$PY" code/test_certificate_contracts.py
"$PY" code/test_magnitude_predicates.py

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
