#!/usr/bin/env python3
"""Deterministic negative test for a winding-two octagon realization.

This is not a proof about all winding-two chambers.  It records that the simplest
exact realization of the combinatorial obstruction is not an NJC counterexample:
its Jacobian determinant becomes negative on a rational grid point.
"""
from __future__ import annotations
import json,math
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
A=[[Q(2),Q(-4),Q(5),Q(-1)],[Q(3),Q(-4),Q(-3),Q(2)]]
B=[[Q(0),Q(1)],[Q(-1),Q(2)],[Q(-2),Q(1)],[Q(-1),Q(-1)]]
c=[Q(0)]*4

def det2(u,v):return u[0]*v[1]-u[1]*v[0]
def sp(t):
    t=float(t);e=math.exp(-abs(t));return e/(1+e)**2
cols=[(A[0][i],A[1][i]) for i in range(4)]
m={(i,j):det2(cols[i],cols[j])*det2(B[i],B[j]) for i in range(4) for j in range(i+1,4)}
def jac(x,y):
 z=[B[i][0]*x+B[i][1]*y for i in range(4)];r=[sp(t) for t in z]
 return sum(float(v)*r[i]*r[j] for (i,j),v in m.items())
# rational point near the deterministic grid minimum
p=(Q(-11,10),Q(-57,20));val=jac(*p)
out={'A':A,'B':B,'mixed_coefficients':{f'{i}{j}':v for (i,j),v in m.items()},'probe_point':p,'determinant_at_probe':val,'conclusion':'NOT_A_COUNTEREXAMPLE' if val<0 else 'INCONCLUSIVE','verdict':'PASS' if val<0 else 'FAIL'}
path=ROOT/'results/structural/n4_winding_obstruction_probe.json';path.write_text(json.dumps(out,default=lambda o:str(o),indent=2)+'\n')
print(json.dumps(out,default=lambda o:str(o),indent=2));returncode=0 if val<0 else 1
raise SystemExit(returncode)
