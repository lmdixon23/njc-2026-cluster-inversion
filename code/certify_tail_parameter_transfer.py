#!/usr/bin/env python3
"""Rigorous transfer of the fixed-witness tail margin to a parameter cube.

The point certificate established a uniform grid-arc tail margin at R=582 of at
least 0.02732.  This script deliberately downgrades that to the exact rational
M0=0.0273 and proves that every parameter perturbation with common l-infinity
radius delta changes every candidate pair margin by less than M0.

The perturbation bound is analytic and uniform in direction.  It combines:
  * four absolute ridge-rate terms: 4*sqrt(2)*R*delta;
  * four bias absolute values: 4*delta;
  * logarithmic changes of the negative and selected positive minor magnitudes;
  * the arc-Lipschitz allowance change: 4*sqrt(2)*half_arc*delta.
Minor-product changes are bounded from the exact centre matrices by expanding the
2x2 determinants and products.  Transcendentals are checked with mpmath.iv.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv, mp


def load(path):
    raw=Path(path).read_bytes(); d=json.loads(raw.decode(),parse_float=Decimal)
    A=[[F(x) for x in r] for r in d['A']]; B=[[F(x) for x in r] for r in d['B']]
    payload=json.dumps({'A':[[str(x) for x in r] for r in d['A']], 'B':[[str(x) for x in r] for r in d['B']], 'c':[str(x) for x in d['c']]},sort_keys=True,separators=(',',':')).encode()
    return A,B,hashlib.sha256(raw).hexdigest(),hashlib.sha256(payload).hexdigest()

def qiv(q):
    x=mp.mpf(q.numerator)/q.denominator
    e=abs(x)*mp.mpf(2)**(1-mp.prec)*8+mp.mpf(2)**(4-mp.prec)
    return iv.mpf([x-e,x+e])

def detA(A,i,j): return A[0][i]*A[1][j]-A[0][j]*A[1][i]
def detB(B,i,j): return B[i][0]*B[j][1]-B[i][1]*B[j][0]

def det_error_A(A,i,j,d):
    C=abs(A[1][j])+abs(A[0][i])+abs(A[1][i])+abs(A[0][j])
    return C*d+2*d*d

def det_error_B(B,i,j,d):
    C=abs(B[j][1])+abs(B[i][0])+abs(B[i][1])+abs(B[j][0])
    return C*d+2*d*d

def minor_error(A,B,i,j,d):
    da=abs(detA(A,i,j)); db=abs(detB(B,i,j))
    ea=det_error_A(A,i,j,d); eb=det_error_B(B,i,j,d)
    return db*ea+da*eb+ea*eb

def main():
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]
    ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--delta',default='5e-6'); ap.add_argument('--delta-a',default=''); ap.add_argument('--delta-b',default=''); ap.add_argument('--delta-c',default=''); ap.add_argument('--R',type=int,default=582)
    ap.add_argument('--grid',type=int,default=200000); ap.add_argument('--dps',type=int,default=60)
    ap.add_argument('--out',default=''); args=ap.parse_args(); mp.dps=args.dps
    d=F(Decimal(args.delta)); dA=F(Decimal(args.delta_a)) if args.delta_a else d; dB=F(Decimal(args.delta_b)) if args.delta_b else d; dC=F(Decimal(args.delta_c)) if args.delta_c else d; A,B,sha,payload_sha=load(args.witness)
    pairs=list(itertools.combinations(range(4),2)); neg=(2,3)
    m={p:detA(A,*p)*detB(B,*p) for p in pairs}
    e={}
    for p in pairs:
        i,j=p; da=abs(detA(A,i,j)); db=abs(detB(B,i,j)); ea=det_error_A(A,i,j,dA); eb=det_error_B(B,i,j,dB); e[p]=db*ea+da*eb+ea*eb
    sign_ok=all(e[p] < abs(m[p]) for p in pairs)
    # interval upper bound on -log(1-e/|m|)
    omega={}
    for p in pairs:
        ratio=qiv(e[p])/qiv(abs(m[p]))
        omega[p]=-iv.log(1-ratio)
    sqrt2=iv.sqrt(2); half_arc=iv.pi/args.grid*(1+mp.mpf(2)**(-30))
    base=qiv(F(273,10000))
    common=(4*sqrt2*args.R*qiv(dB) + 4*qiv(dC) + 4*sqrt2*half_arc*qiv(dB))
    degrad={p: common+omega[neg]+omega[p] for p in pairs if p!=neg}
    worst_pair=max(degrad,key=lambda p: mp.mpf(degrad[p].b))
    worst=degrad[worst_pair]
    residual=base-worst
    passed=sign_ok and residual.a>0
    out={
      'verdict':'PASS' if passed else 'FAIL','uniform_parameter_radius_inf':args.delta,'delta_A':str(dA),'delta_B':str(dB),'delta_c':str(dC),
      'R':args.R,'grid':args.grid,'base_tail_margin_lower':'0.0273',
      'witness_sha256':sha,'ABC_payload_sha256':payload_sha,
      'minor_sign_preservation':sign_ok,
      'minor_error_bounds':{str(p):str(e[p]) for p in pairs},
      'worst_pair':list(worst_pair),
      'tail_margin_degradation_upper':str(worst.b),
      'tail_margin_residual_lower':str(residual.a),
      'method':'analytic uniform transfer from certified point-tail grid margin'
    }
    text=json.dumps(out,indent=2)+'\n'; print(text,end='')
    if args.out: Path(args.out).write_text(text, newline="\n")
    return 0 if passed else 1
if __name__=='__main__': raise SystemExit(main())
