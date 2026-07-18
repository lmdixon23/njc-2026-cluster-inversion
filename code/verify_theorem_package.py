#!/usr/bin/env python3
"""Fast semantic release verification for the theorem and certificate package."""
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path
from report_validation import (
    ReportValidationError,
    file_sha256,
    load_json,
    validate_combined,
    validate_point_tail,
    validate_point_tail_arb,
)
ROOT=Path(__file__).resolve().parents[1]

def run(script,*args):
    cmd=[sys.executable,str(ROOT/'code'/script),*map(str,args)]; print('+',' '.join(cmd),flush=True)
    subprocess.run(cmd,cwd=ROOT,check=True)

def main():
    run('verify_witness_polygon.py','--json-out','results/verification_report.json')
    run('verify_sign_only_counterexample.py')
    run('verify_comparison_examples.py')
    try:
        validate_combined(ROOT/'results'/'witness_common_5e-6', ROOT/'config'/'witness_common.json')
        validate_combined(ROOT/'results'/'witness_anisotropic', ROOT/'config'/'witness_anisotropic.json')
        witness=ROOT/'data'/'witness.json'
        point=ROOT/'results'/'point_tail_recomputed.json'
        validate_point_tail(load_json(point), file_sha256(witness), 582, 200000)
        validate_point_tail_arb(ROOT/'results'/'point_tail_arb.json', point, witness, 582, 200000)
    except ReportValidationError as exc:
        raise SystemExit(f'theorem package semantic validation failed: {exc}') from exc
    print('[PASS] semantically validated common radius, anisotropic box, and complete Arb replays')
    return 0
if __name__=='__main__': raise SystemExit(main())
