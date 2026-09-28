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
  * optional positive grid slack: 4*sqrt(2)*half_arc*delta.
The retained base margin is already uniform in angle, so the transfer does not
need to recompute an arc-Lipschitz allowance. Minor-product changes are bounded from the exact centre matrices by expanding the
2x2 determinants and products.  Transcendentals are checked with mpmath.iv.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv, mp
from report_validation import (BASE_TAIL_MARGIN, require, validate_planar_four, validated_tail_base)


def load(path):
    raw=Path(path).read_bytes(); d=json.loads(raw.decode(),parse_float=Decimal)
    A=[[F(x) for x in r] for r in d['A']]; B=[[F(x) for x in r] for r in d['B']]
    c=[F(x) for x in d['c']]
    validate_planar_four(A, B, c)
    payload=json.dumps({'A':[[str(x) for x in r] for r in d['A']], 'B':[[str(x) for x in r] for r in d['B']], 'c':[str(x) for x in d['c']]},sort_keys=True,separators=(',',':')).encode()
    return A,B,c,hashlib.sha256(raw).hexdigest(),hashlib.sha256(payload).hexdigest()

def qiv(q):
    return iv.mpf(q.numerator)/iv.mpf(q.denominator)

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


def worst_positive_pair(minors, errors, negative):
    # -log(1-r) is increasing for 0<=r<1. Select its largest exact
    # rational argument, not a rounded ranking of interval endpoints.
    return max((p for p in minors if p != negative), key=lambda p: errors[p]/abs(minors[p]))

def main():
    ap=argparse.ArgumentParser(); root=Path(__file__).resolve().parents[1]
    ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--delta',default='5e-6'); ap.add_argument('--delta-a',default=''); ap.add_argument('--delta-b',default=''); ap.add_argument('--delta-c',default=''); ap.add_argument('--R',type=int,default=582)
    ap.add_argument('--grid',type=int,default=200000); ap.add_argument('--dps',type=int,default=60)
    ap.add_argument('--out',default='')
    ap.add_argument('--point-report', default=''); ap.add_argument('--point-arb-report', default='')
    args=ap.parse_args()
    require(args.dps >= 30, 'precision must be at least 30 decimal digits')
    mp.dps=iv.dps=args.dps
    d=F(Decimal(args.delta)); dA=F(Decimal(args.delta_a)) if args.delta_a else d; dB=F(Decimal(args.delta_b)) if args.delta_b else d; dC=F(Decimal(args.delta_c)) if args.delta_c else d; A,B,c,sha,payload_sha=load(args.witness)
    require(all(radius >= 0 for radius in (d, dA, dB, dC)), 'parameter radii must be nonnegative')
    base_binding=validated_tail_base(Path(args.witness),args.R,args.grid,args.point_report,args.point_arb_report)
    pairs=list(itertools.combinations(range(4),2)); neg=(2,3)
    m={p:detA(A,*p)*detB(B,*p) for p in pairs}
    require([p for p in pairs if m[p]<0] == [neg] and sum(m[p]>0 for p in pairs)==5,
            'expected five positive minors and negative pair (2,3)')
    e={}
    for p in pairs:
        i,j=p; da=abs(detA(A,i,j)); db=abs(detB(B,i,j)); ea=det_error_A(A,i,j,dA); eb=det_error_B(B,i,j,dB); e[p]=db*ea+da*eb+ea*eb
    sign_ok=all(e[p] < abs(m[p]) for p in pairs)
    require(sign_ok, 'minor signs are not uniformly preserved')
    # A positive R*g-K transfers to all larger radii only with g>0. A
    # uniform K>0 is sufficient here and is checked for the whole box.
    for p in pairs:
        if p == neg: continue
        k_lower=iv.log(16*qiv(abs(m[neg])-e[neg])/qiv(m[p]+e[p]))
        k_lower+=sum(qiv(max(abs(c[i])-dC,F(0))) for i in (*neg,*p))
        require(bool(k_lower>0), 'cannot establish positive tail K throughout parameter box')
    # interval upper bound on -log(1-e/|m|)
    omega={}
    for p in pairs:
        ratio=qiv(e[p])/qiv(abs(m[p]))
        omega[p]=-iv.log(1-ratio)
    sqrt2=iv.sqrt(2); half_arc=iv.pi/args.grid*(1+mp.mpf(2)**(-30))
    base=qiv(BASE_TAIL_MARGIN)
    # The base is already uniform over each angular arc. Directionwise
    # transfer needs R*Delta(g)+Delta(K); the extra positive arc term is
    # conservative slack, not a recomputation of grid Lipschitz allowances.
    common=(4*sqrt2*args.R*qiv(dB) + 4*qiv(dC) + 4*sqrt2*half_arc*qiv(dB))
    degrad={p: common+omega[neg]+omega[p] for p in pairs if p!=neg}
    worst_pair=worst_positive_pair(m,e,neg)
    worst=degrad[worst_pair]
    residual=base-worst
    passed=sign_ok and residual.a>0
    out={
      'schema':'tail-transfer-v2','base_evidence':base_binding,
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
