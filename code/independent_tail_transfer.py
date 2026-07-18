#!/usr/bin/env python3
"""Independent Decimal verification of the analytic tail-margin transfer.

This script re-derives minor perturbation bounds using Decimal arithmetic and
checks the explicit degradation formula without importing the primary helpers.
"""
from __future__ import annotations
import argparse,itertools,json
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

P=80

def D(x): return x if isinstance(x,Decimal) else Decimal(str(x))
def up(fn):
    with localcontext() as c: c.prec=P; c.rounding=ROUND_CEILING; v=fn()
    return v + D('1e-70')*(1+abs(v))
def down(fn):
    with localcontext() as c: c.prec=P; c.rounding=ROUND_FLOOR; v=fn()
    return v - D('1e-70')*(1+abs(v))
def sqrt_up(x): return up(lambda: x.sqrt())
def ln_up(x): return up(lambda: x.ln())

def detA(A,i,j): return A[0][i]*A[1][j]-A[0][j]*A[1][i]
def detB(B,i,j): return B[i][0]*B[j][1]-B[i][1]*B[j][0]
def det_err(M,i,j,d):
    C=abs(M[1][j])+abs(M[0][i])+abs(M[1][i])+abs(M[0][j])
    return up(lambda: C*d+2*d*d)
def det_err_rows(M,i,j,d):
    C=abs(M[j][1])+abs(M[i][0])+abs(M[i][1])+abs(M[j][0])
    return up(lambda: C*d+2*d*d)

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(); ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--delta',default='0.000005'); ap.add_argument('--delta-a',default=''); ap.add_argument('--delta-b',default=''); ap.add_argument('--delta-c',default='')
    ap.add_argument('--R',default='582'); ap.add_argument('--grid',default='200000'); ap.add_argument('--base-margin',default='0.0273'); ap.add_argument('--json-out',default='')
    a=ap.parse_args(); d=D(a.delta); da=D(a.delta_a) if a.delta_a else d; db=D(a.delta_b) if a.delta_b else d; dc=D(a.delta_c) if a.delta_c else d
    data=json.loads(Path(a.witness).read_text(),parse_float=Decimal,parse_int=Decimal); A=data['A']; B=data['B']
    pairs=list(itertools.combinations(range(len(B)),2)); neg=(2,3)
    m={p:detA(A,*p)*detB(B,*p) for p in pairs}; err={}
    for p in pairs:
        i,j=p; ma=abs(detA(A,i,j)); mb=abs(detB(B,i,j)); ea=det_err(A,i,j,da); eb=det_err_rows(B,i,j,db)
        err[p]=up(lambda mb=mb,ma=ma,ea=ea,eb=eb: mb*ea+ma*eb+ea*eb)
    signs=all(err[p] < abs(m[p]) for p in pairs)
    omega={p: -down(lambda p=p: (D(1)-err[p]/abs(m[p])).ln()) for p in pairs}
    sqrt2=sqrt_up(D(2)); pi_up=D('3.1415926535897932384626433832795028841971693993751058209749445923078164063')
    half_arc=up(lambda: pi_up/D(a.grid)*(D(1)+D(2)**D(-30)))
    common=up(lambda: 4*sqrt2*D(a.R)*db + 4*dc + 4*sqrt2*half_arc*db)
    degradation={p:up(lambda p=p: common+omega[neg]+omega[p]) for p in pairs if p!=neg}
    worst=max(degradation,key=degradation.get); residual=down(lambda: D(a.base_margin)-degradation[worst])
    out={'verdict':'PASS' if signs and residual>0 else 'FAIL','delta_A':str(da),'delta_B':str(db),'delta_c':str(dc),'minor_sign_preservation':signs,'minor_error_upper':{str(k):str(v) for k,v in err.items()},'log_loss_upper':{str(k):str(v) for k,v in omega.items()},'common_rate_bias_arc_loss_upper':str(common),'worst_pair':list(worst),'degradation_upper':str(degradation[worst]),'residual_lower':str(residual),'arithmetic':'independent Decimal, 80 digits, outward padding'}
    print(json.dumps(out,indent=2))
    if a.json_out: Path(a.json_out).write_text(json.dumps(out,indent=2)+'\n')
    raise SystemExit(0 if out['verdict']=='PASS' else 1)
if __name__=='__main__': main()
