#!/usr/bin/env python3
"""Diagnostic replay of exported leaf paths with custom Decimal intervals.

The primary certificate exports a prefix-free set of binary split paths. This checker
reconstructs each dyadic leaf from the root, verifies exact partition volume and
prefix-freeness, then diagnostically reevaluates every leaf with a separately
written Decimal interval implementation. The primary interval certificate is
the rigorous source of the compact positivity claim.
"""
from __future__ import annotations
import argparse,json,math,multiprocessing as mp,time
from decimal import Decimal
from pathlib import Path
from independent_parameter_checker import load,params,minors,det_interval
from decimal_interval import D

_DATA=None; _M=None; _B=None; _C=None; _R=None



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

def init_worker(witness,da,db,dc,R):
    global _DATA,_M,_B,_C,_R
    _DATA=load(Path(witness)); A,_B,_C=params(_DATA,D(da),D(db),D(dc)); _M=minors(A,_B); _R=Decimal(str(R))

def path_box(path:str):
    x1=y1=-_R; x2=y2=_R
    for bit in path:
        w=x2-x1; h=y2-y1
        if w>=h:
            m=(x1+x2)/2
            if bit=='0': x2=m
            else: x1=m
        else:
            m=(y1+y2)/2
            if bit=='0': y2=m
            else: y1=m
    return x1,y1,x2,y2

def check_one(path):
    box=path_box(path); iv=det_interval(box,_B,_C,_M)
    return path, str(iv.lo), str(iv.hi), iv.lo>0

def prefix_free(paths):
    s=sorted(paths)
    return all(not s[i+1].startswith(s[i]) for i in range(len(s)-1))

def exact_volume(paths,R):
    # Each binary split halves area, so normalized volume is sum 2^-depth.
    # Use integer common denominator at max depth.
    md=max(map(len,paths)); num=sum(1<<(md-len(p)) for p in paths)
    return num==(1<<md), {'max_depth':md,'normalized_numerator':str(num),'normalized_denominator':str(1<<md),'root_area':str((Decimal(2*R))**2)}

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--witness',default=str(root/'data'/'witness.json'))
    ap.add_argument('--paths',default=str(root/'verification'/'compact_leaf_paths.json'))
    ap.add_argument('--delta',default='0.000005'); ap.add_argument('--delta-a',default=''); ap.add_argument('--delta-b',default=''); ap.add_argument('--delta-c',default='')
    ap.add_argument('--workers',type=int,default=max(1,min(8,mp.cpu_count())))
    ap.add_argument('--chunksize',type=int,default=256); ap.add_argument('--json-out',default='')
    a=ap.parse_args(); d=a.delta; da=a.delta_a or d; db=a.delta_b or d; dc=a.delta_c or d
    rec=json.load(open(a.paths)); paths=rec['paths']; R=int(rec['R'])
    pf=prefix_free(paths); vol_ok,vol=exact_volume(paths,R)
    t=time.time(); bad=[]; minlo=None; checked=0
    with mp.Pool(a.workers,initializer=init_worker,initargs=(a.witness,da,db,dc,R)) as pool:
        for path,lo,hi,ok in pool.imap_unordered(check_one,paths,chunksize=a.chunksize):
            checked+=1; dlo=Decimal(lo); minlo=dlo if minlo is None else min(minlo,dlo)
            if not ok and len(bad)<10: bad.append({'path':path,'interval':[lo,hi]})
    out={'verdict':'PASS' if pf and vol_ok and not bad and checked==len(paths) else 'FAIL','arithmetic':'custom Decimal interval diagnostic','delta_A':da,'delta_B':db,'delta_c':dc,'leaf_count':len(paths),'checked':checked,'prefix_free':pf,'exact_volume_closure':vol_ok,'volume':vol,'minimum_leaf_lower':str(minlo),'failures':bad,'elapsed_seconds':time.time()-t,'workers':a.workers}
    print(json.dumps(out,indent=2))
    if a.json_out: Path(a.json_out).write_text(json.dumps(_strip_timing(out),indent=2)+'\n', newline="\n")
    raise SystemExit(0 if out['verdict']=='PASS' else 1)
if __name__=='__main__': main()
