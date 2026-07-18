#!/usr/bin/env python3
"""Independent Decimal-interval replay of the compact parameter certificate.

This checker shares no arithmetic helper with certify_parameter_radius.py.  It
reads the same witness and common/anisotropic parameter box, recomputes all six
minor-product intervals, and replays an exhaustive fail-closed branch-and-bound
on the compact square. It is intentionally slower than the primary certificate.
"""
from __future__ import annotations
import argparse, itertools, json, time
from decimal import Decimal
from pathlib import Path
from decimal_interval import DI, D, sigma_prime


def load(path: Path):
    return json.loads(path.read_text(encoding='utf-8'), parse_float=Decimal, parse_int=Decimal)


def box(v: Decimal, d: Decimal) -> DI:
    return DI(v-d, v+d)


def params(data, da, db, dc):
    A=[[box(v,da) for v in row] for row in data['A']]
    B=[[box(v,db) for v in row] for row in data['B']]
    c=[box(v,dc) for v in data['c']]
    return A,B,c


def minors(A,B):
    out={}
    n=len(B)
    for i,j in itertools.combinations(range(n),2):
        dA=A[0][i]*A[1][j]-A[0][j]*A[1][i]
        dB=B[i][0]*B[j][1]-B[i][1]*B[j][0]
        out[(i,j)]=dA*dB
    return out


def det_interval(rect,B,c,M):
    x1,y1,x2,y2=rect
    X=DI(x1,x2); Y=DI(y1,y2)
    sp=[]
    for i in range(len(B)):
        sp.append(sigma_prime(B[i][0]*X+B[i][1]*Y+c[i]))
    acc=DI.point(0)
    for (i,j),m in M.items():
        acc=acc+m*sp[i]*sp[j]
    return acc


def certify(data, da, db, dc, R, hmin, budget, progress):
    A,B,c=params(data,da,db,dc)
    M=minors(A,B)
    signs={str(k): ('+' if v.lo>0 else '-' if v.hi<0 else '?') for k,v in M.items()}
    if '?' in signs.values():
        return {'verdict':'FAIL','failure':'minor ambiguity','minor_signs':signs}
    root=(Decimal(-R),Decimal(-R),Decimal(R),Decimal(R))
    stack=[root]; processed=leaves=0; failures=[]; minlo=None; t0=time.time()
    while stack and processed<budget:
        rect=stack.pop(); processed+=1
        iv=det_interval(rect,B,c,M)
        minlo=iv.lo if minlo is None else min(minlo,iv.lo)
        if iv.lo>0:
            leaves+=1
        else:
            x1,y1,x2,y2=rect; w=x2-x1; h=y2-y1
            if w<hmin and h<hmin:
                failures.append({'rect':[str(z) for z in rect],'det':[str(iv.lo),str(iv.hi)]})
                if len(failures)>=3: break
            elif w>=h:
                m=(x1+x2)/2; stack.append((x1,y1,m,y2)); stack.append((m,y1,x2,y2))
            else:
                m=(y1+y2)/2; stack.append((x1,y1,x2,m)); stack.append((x1,m,x2,y2))
        if progress and processed%progress==0:
            print(f'DECIMAL replay processed={processed:,} stack={len(stack):,}',flush=True)
    exhausted=bool(stack)
    return {
      'verdict':'PASS' if not failures and not exhausted else 'FAIL',
      'arithmetic':'custom Decimal intervals, precision=70, explicit endpoint padding',
      'delta_A':str(da),'delta_B':str(db),'delta_c':str(dc),
      'R':R,'hmin':str(hmin),'processed':processed,'certified_leaves':leaves,
      'budget_exhausted':exhausted,'unresolved':failures,'smallest_seen_lower':str(minlo),
      'minor_signs':signs,'elapsed_seconds':time.time()-t0,
    }


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--delta',default='0.000005')
    ap.add_argument('--delta-a',default='')
    ap.add_argument('--delta-b',default='')
    ap.add_argument('--delta-c',default='')
    ap.add_argument('--R',type=int,default=582)
    ap.add_argument('--hmin',default='0.008')
    ap.add_argument('--budget',type=int,default=1000000)
    ap.add_argument('--progress',type=int,default=50000)
    ap.add_argument('--json-out',default='')
    a=ap.parse_args(); d=D(a.delta)
    da=D(a.delta_a) if a.delta_a else d; db=D(a.delta_b) if a.delta_b else d; dc=D(a.delta_c) if a.delta_c else d
    res=certify(load(Path(a.witness)),da,db,dc,a.R,D(a.hmin),a.budget,a.progress)
    print(json.dumps(res,indent=2))
    if a.json_out: Path(a.json_out).write_text(json.dumps(res,indent=2)+'\n')
    raise SystemExit(0 if res['verdict']=='PASS' else 1)

if __name__=='__main__': main()
