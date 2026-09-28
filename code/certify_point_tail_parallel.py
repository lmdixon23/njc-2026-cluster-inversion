#!/usr/bin/env python3
"""Parallel rigorous recomputation of the point-witness fixed-radius tail margin.

This is a provenance-clean replacement for stale external-session reports.  It
reads the canonical witness, builds the same validated interval constants as
``certify_parameter_radius.py`` at delta=0, partitions the angular grid into
independent contiguous chunks, and returns the global minimum of the per-arc
best domination margin.  Each worker uses mpmath.iv outward interval arithmetic.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, multiprocessing as mpool, time
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv, mp
from report_validation import validate_planar_four, validate_fixed_context, require




def _strip_timing(payload):
    """Return a copy with wall-clock fields removed, recursively.

    Persisted artifacts are content-addressed in SHA256SUMS.txt, so their bytes
    must depend only on the mathematics.  Timings stay on stdout.
    """
    if isinstance(payload, dict):
        return {k: _strip_timing(v) for k, v in payload.items()
                if k != "elapsed_seconds"}
    if isinstance(payload, list):
        return [_strip_timing(v) for v in payload]
    return payload

def load(path: Path):
    raw=path.read_bytes(); d=json.loads(raw.decode('utf-8'),parse_float=Decimal)
    A=[[F(x) for x in r] for r in d['A']]
    B=[[F(x) for x in r] for r in d['B']]
    c=[F(x) for x in d['c']]
    validate_planar_four(A, B, c)
    return A,B,c,hashlib.sha256(raw).hexdigest()

def qiv(q:F):
    x=mp.mpf(q.numerator)/mp.mpf(q.denominator)
    e=abs(x)*mp.mpf(2)**(1-mp.prec)*8+mp.mpf(2)**(4-mp.prec)
    return iv.mpf([x-e,x+e])

def iabs(X):
    if X.a>=0:return X
    if X.b<=0:return -X
    return iv.mpf([0,max(-X.a,X.b)])

def norm_up(v):
    x=mp.mpf(v[0].numerator)/v[0].denominator
    y=mp.mpf(v[1].numerator)/v[1].denominator
    return mp.sqrt(x*x+y*y)*(1+mp.mpf(2)**(-20))

def minors(A,B):
    validate_planar_four(A, B)
    out={}
    for i,j in itertools.combinations(range(4),2):
        da=A[0][i]*A[1][j]-A[0][j]*A[1][i]
        db=B[i][0]*B[j][1]-B[i][1]*B[j][0]
        out[(i,j)]=da*db
    return out

def constants(path,dps):
    mp.dps=dps
    A,B,c,sha=load(Path(path)); m=minors(A,B)
    pos=[p for p,v in m.items() if v>0]; neg=[p for p,v in m.items() if v<0]
    if neg!=[(2,3)]: raise RuntimeError(f'unexpected negative set {neg}')
    I=neg[0]; mn=-m[I]
    K={}
    L={}
    for J in pos:
        K[J]=iv.log(16*qiv(mn)/qiv(m[J]))+sum(qiv(abs(c[k])) for k in (I[0],I[1],J[0],J[1]))
        L[J]=sum(norm_up(B[k]) for k in (I[0],I[1],J[0],J[1]))
    return B,m,pos,I,K,L,sha

def worker(args):
    path,R,Mgrid,start,end,dps=args
    require(type(R) is int and R > 0 and type(Mgrid) is int and Mgrid > 0
            and 0 <= start < end <= Mgrid and dps >= 30, 'invalid angular worker context')
    B,m,pos,I,K,L,sha=constants(path,dps)
    half_arc=mp.pi/Mgrid*(1+mp.mpf(2)**(-30))
    dth=2*iv.pi/Mgrid
    worst=mp.inf; worstk=-1; worstpair=None
    def rho(i,ct,st): return iabs(qiv(B[i][0])*ct+qiv(B[i][1])*st)
    for k in range(start,end):
        th=dth*k; ct,st=iv.cos(th),iv.sin(th)
        rn=rho(I[0],ct,st)+rho(I[1],ct,st)
        best=-mp.inf; pair=None
        for J in pos:
            g=(rn-rho(J[0],ct,st)-rho(J[1],ct,st))-iv.mpf(L[J]*half_arc)
            if mp.mpf(g.a)<=0: continue
            margin=mp.mpf((R*g-K[J]).a)
            if margin>best: best=margin; pair=J
        if best<worst: worst=best; worstk=k; worstpair=pair
        if not best>0:
            return {'passed':False,'worst':str(worst),'worstk':worstk,'worst_pair':worstpair,'start':start,'end':end}
    return {'passed':True,'worst':str(worst),'worstk':worstk,'worst_pair':worstpair,'start':start,'end':end}

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(); ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--R',type=int,default=582); ap.add_argument('--grid',type=int,default=200000)
    ap.add_argument('--workers',type=int,default=max(1,min(8,mpool.cpu_count()))); ap.add_argument('--dps',type=int,default=40)
    ap.add_argument('--json-out',default=str(root/'results'/'point_tail_recomputed.json'))
    a=ap.parse_args(); A,B,c,sha=load(Path(a.witness)); m=minors(A,B)
    validate_fixed_context(sha, a.R, a.grid)
    require(a.workers > 0 and a.dps >= 30, 'invalid workers or precision')
    chunks=[]
    for w in range(a.workers):
        lo=a.grid*w//a.workers; hi=a.grid*(w+1)//a.workers
        if lo<hi: chunks.append((a.witness,a.R,a.grid,lo,hi,a.dps))
    t=time.time()
    ctx=mpool.get_context('spawn')
    with ctx.Pool(len(chunks)) as pool: parts=pool.map(worker,chunks)
    passed=all(p['passed'] for p in parts)
    worstpart=min(parts,key=lambda p: mp.mpf(p['worst']))
    out={
      'delta':'0','delta_A':'0','delta_B':'0','delta_c':'0','witness_sha256':sha,
      'R':a.R,'grid':a.grid,'minor_signs':{'positive':[list(p) for p,v in m.items() if v>0],'negative':[list(p) for p,v in m.items() if v<0],'ambiguous':[]},
      'tail':{'passed':passed,'worst_margin':worstpart['worst'],'worst_grid_index':worstpart['worstk'],'worst_pair':worstpart['worst_pair'],'workers':len(chunks),'chunks':parts,'elapsed_seconds':time.time()-t},
      'verdict':'PASS' if passed and mp.mpf(worstpart['worst'])>0 else 'FAIL',
      'method':'parallel mpmath.iv fixed-radius angular-grid certificate; exact rational witness input'
    }
    text=json.dumps(out,indent=2)+'\n'; print(text,end=''); Path(a.json_out).write_text(json.dumps(_strip_timing(out),indent=2)+'\n', newline="\n")
    return 0 if out['verdict']=='PASS' else 1
if __name__=='__main__': raise SystemExit(main())
