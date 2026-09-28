#!/usr/bin/env python3
"""Rigorous scalar certificate for an analytic logistic-Jacobian dominance theorem.

The theorem uses one negative mixed coefficient m_23 and the two positive terms
m_02, m_03.  If b0 = alpha*b2 + beta*b3 with |alpha|+|beta|<1, and weights
mu+lambda=1 satisfy mu>=|alpha|, lambda>=|beta|, then convexity of
psi(t)=2 log cosh(t/2) gives a global lower bound on
m02*rho0/rho3 + m03*rho0/rho2.

This script derives every algebraic quantity exactly from decimal JSON input and
uses mpmath.iv only for the final transcendental scalar inequality.
"""
from __future__ import annotations
import argparse, hashlib, json
from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
import mpmath as mp
from report_validation import validate_planar_four


def det(u,v): return u[0]*v[1]-u[1]*v[0]
def qs(x: Q) -> str: return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)
def ivq(x: Q): return mp.iv.mpf([str(x.numerator),str(x.numerator)])/mp.iv.mpf([str(x.denominator),str(x.denominator)])
def iv_bounds(x): return str(x.a), str(x.b)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--witness',default='data/witness.json')
    ap.add_argument('--mu',default='51/1000',help='weight on m03*rho0/rho2')
    ap.add_argument('--json-out',default='results/structural/analytic_determinant.json')
    ns=ap.parse_args()
    raw=Path(ns.witness).read_bytes()
    obj=json.loads(raw.decode(),parse_float=Decimal)
    A=[[Q(x) for x in row] for row in obj['A']]
    B=[[Q(x) for x in row] for row in obj['B']]
    c=[Q(x) for x in obj['c']]
    validate_planar_four(A, B, c)
    cols=[(A[0][i],A[1][i]) for i in range(len(B))]
    m={(i,j):det(cols[i],cols[j])*det(B[i],B[j]) for i in range(4) for j in range(i+1,4)}
    D=det(B[2],B[3])
    if D == 0: raise SystemExit('rows b2,b3 are dependent')
    alpha=det(B[0],B[3])/D
    beta=det(B[2],B[0])/D
    delta=c[0]-alpha*c[2]-beta*c[3]
    a,b=abs(alpha),abs(beta)
    r=Q(1)-a-b
    mu=Q(ns.mu); lam=Q(1)-mu
    hypotheses={
      'r_positive':r>0,'mu_positive':mu>0,'lambda_positive':lam>0,
      'mu_at_least_abs_alpha':mu>=a,'lambda_at_least_abs_beta':lam>=b,
      'm01_nonnegative':m[0,1]>=0,'m02_positive':m[0,2]>0,
      'm03_positive':m[0,3]>0,'m12_nonnegative':m[1,2]>=0,
      'm13_nonnegative':m[1,3]>=0,'m23_negative':m[2,3]<0,
    }
    algebraic=all(hypotheses.values())
    mp.iv.dps=90
    R, MU, LA = ivq(r),ivq(mu),ivq(lam)
    DEL=ivq(delta)
    M02,M03,M23=ivq(m[0,2]),ivq(m[0,3]),ivq(-m[2,3])
    t=abs(DEL/R)
    psi=2*mp.iv.log((mp.iv.exp(t/2)+mp.iv.exp(-t/2))/2)
    log_lower = LA*(mp.iv.log(M02)-mp.iv.log(LA)) + MU*(mp.iv.log(M03)-mp.iv.log(MU)) - R*psi
    lower=mp.iv.exp(log_lower)
    margin=lower-M23
    passed=bool(algebraic and margin>0)
    report={
      'witness_sha256':hashlib.sha256(raw).hexdigest(),
      'weights':{'mu':qs(mu),'lambda':qs(lam)},
      'row_decomposition':{'alpha':qs(alpha),'beta':qs(beta),'delta':qs(delta),'r':qs(r)},
      'mixed_coefficients':{f'{i}{j}':qs(v) for (i,j),v in m.items()},
      'hypotheses':hypotheses,
      'algebraic_hypotheses':algebraic,
      'log_lower_interval':iv_bounds(log_lower),
      'dominance_lower_interval':iv_bounds(lower),
      'negative_coefficient_interval':iv_bounds(M23),
      'margin_interval':iv_bounds(margin),
      'verdict':'PASS' if passed else 'FAIL'
    }
    out=Path(ns.json_out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n', newline="\n")
    print(json.dumps(report,indent=2))
    return 0 if passed else 1
if __name__=='__main__': raise SystemExit(main())
