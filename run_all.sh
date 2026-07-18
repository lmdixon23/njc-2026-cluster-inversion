#!/usr/bin/env bash
# Reproduce the maintained deterministic claims in mathematical order.
set -euo pipefail
export PYTHONHASHSEED=0

PY=python3
"$PY" -c "" 2>/dev/null || PY=python
"$PY" -c "" 2>/dev/null || { echo "no working Python interpreter found"; exit 1; }

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

printf '%s\n' '== Independent analytic determinant reconstruction =='
"$PY" code/check_analytic_det.py

printf '%s\n' '== Independent witness reconstruction =='
"$PY" code/check_witness_core.py

printf '%s\n' '== Independent exact polygon reconstruction =='
"$PY" code/polygon_indep.py

printf '%s\n' '== Final canonical-file integrity =='
"$PY" code/verify_sha256_manifest.py

printf '%s\n' 'VERDICT: ALL MAINTAINED CHECKS PASSED'
