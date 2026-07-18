#!/usr/bin/env python3
"""Fail-closed runner for the structural theorem checks."""
from __future__ import annotations
import json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def run(name,*args):
    cmd=[sys.executable,str(ROOT/'code'/name),*args]
    print('+',' '.join(cmd),flush=True);subprocess.run(cmd,cwd=ROOT,check=True)

def req(path):
    p=ROOT/path;obj=json.load(open(p));
    if obj.get('verdict')!='PASS':raise SystemExit(f'{path}: not PASS')
    print('[PASS]',path)

def main():
    run('certify_analytic_determinant.py')
    run('verify_analytic_family_example.py')
    run('polygon_winding.py')
    run('magnitude_chambers.py')
    run('verify_magnitude_chamber_refinement.py')
    run('probe_n4_winding_obstruction.py')
    for p in ['results/structural/analytic_determinant.json','results/structural/analytic_family_example.json','results/structural/winding_spectra.json','results/structural/magnitude_chambers.json','results/structural/magnitude_chamber_refinement_witness.json','results/structural/n4_winding_obstruction_probe.json']:
        req(p)
    print('[PASS] structural theorem checks')
    return 0
if __name__=='__main__':raise SystemExit(main())
