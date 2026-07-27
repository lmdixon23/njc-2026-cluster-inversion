#!/usr/bin/env python3
"""Exact geometry plus analytic determinant certificate for a non-witness family member."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from asymptotic_polygon import load_exact_json, build_polygon, polygon_simplicity, strict_convexity_by_turns, mixed_minors, best_integer_functional, monotone_one_sided, jsonable

def main():
    path=ROOT/'data/analytic_family_example.json'
    d=load_exact_json(path); poly=build_polygon(d['A'],d['B'])
    simple,fail=polygon_simplicity(poly.vertices); convex,_=strict_convexity_by_turns(poly.edges)
    q=best_integer_functional(poly.half_edges,limit=100)
    one=monotone_one_sided(poly.prefixes,q) if q else {'passed':False}
    minors=mixed_minors(d['A'],d['B'])
    cert=ROOT/'results/structural/analytic_family_example_determinant.json'
    if cert.exists(): cert.unlink()
    proc=subprocess.run([sys.executable,str(ROOT/'code/certify_analytic_determinant.py'),'--witness',str(path),'--mu','1/10','--json-out',str(cert)],cwd=ROOT,check=True)
    det_report=json.load(open(cert))
    out={'polygon_simple':simple,'strictly_convex':convex,'failures':fail,'q':q,'one_sided':one,'mixed_minors':minors,'determinant_certificate':det_report,'verdict':'PASS' if simple and not convex and det_report['verdict']=='PASS' else 'FAIL'}
    op=ROOT/'results/structural/analytic_family_example.json';op.parent.mkdir(parents=True,exist_ok=True);op.write_text(json.dumps(jsonable(out),indent=2)+'\n', newline="\n")
    print(json.dumps(jsonable(out),indent=2));return 0 if out['verdict']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
