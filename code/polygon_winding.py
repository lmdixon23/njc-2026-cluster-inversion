#!/usr/bin/env python3
"""Exact winding numbers and vertical-decomposition spectra for rational polygons."""
from __future__ import annotations
import argparse,json,sys
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from asymptotic_polygon import load_exact_json,build_polygon,det,sub,on_segment,jsonable

def winding_number(vertices,p):
    x,y=p; w=0
    for i,a in enumerate(vertices):
        b=vertices[(i+1)%len(vertices)]
        if on_segment(p,a,b): return None
        cross=det(sub(b,a),sub(p,a))
        if a[1] <= y < b[1] and cross>0: w+=1
        elif b[1] <= y < a[1] and cross<0: w-=1
    return w

def line_intersection(a,b,c,d):
    u=sub(b,a);v=sub(d,c);den=det(u,v)
    if den==0:return None
    t=det(sub(c,a),v)/den
    s=det(sub(c,a),u)/den
    if 0<=t<=1 and 0<=s<=1:return (a[0]+t*u[0],a[1]+t*u[1])
    return None

def spectrum(vertices):
    xs={v[0] for v in vertices}
    n=len(vertices)
    for i in range(n):
      for j in range(i+1,n):
        if j==i+1 or (i==0 and j==n-1):continue
        p=line_intersection(vertices[i],vertices[(i+1)%n],vertices[j],vertices[(j+1)%n])
        if p is not None:xs.add(p[0])
    sx=sorted(xs)
    if not sx:return {}
    span=max(Q(1),sx[-1]-sx[0])
    sample_x=[sx[0]-span,sx[-1]+span]+[(a+b)/2 for a,b in zip(sx,sx[1:])]
    counts={};samples={}
    for x in sample_x:
        ys=[]
        for i,a in enumerate(vertices):
            b=vertices[(i+1)%n]
            lo,hi=sorted((a[0],b[0]))
            if lo < x < hi:
                t=(x-a[0])/(b[0]-a[0]);ys.append(a[1]+t*(b[1]-a[1]))
        ys=sorted(set(ys))
        if ys:
            yspan=max(Q(1),ys[-1]-ys[0]); sy=[ys[0]-yspan,ys[-1]+yspan]+[(a+b)/2 for a,b in zip(ys,ys[1:])]
        else: sy=[Q(0)]
        for y in sy:
            w=winding_number(vertices,(x,y))
            if w is not None:
                counts[w]=counts.get(w,0)+1;samples.setdefault(w,(x,y))
    return {'winding_values':sorted(counts),'sample_counts':counts,'samples':samples,'max_abs_winding':max(map(abs,counts)) if counts else 0}

def obstruction_polygon():
    half=[(Q(2),Q(3)),(Q(-4),Q(-4)),(Q(5),Q(-3)),(Q(-1),Q(2))]
    edges=half+[(-x,-y) for x,y in half];v=[(Q(0),Q(0))]
    for e in edges:v.append((v[-1][0]+e[0],v[-1][1]+e[1]))
    return v[:-1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--witness',default='data/witness.json');ap.add_argument('--json-out',default='results/structural/winding_spectra.json');ns=ap.parse_args()
    d=load_exact_json(ns.witness);wp=build_polygon(d['A'],d['B']).vertices
    ob=obstruction_polygon()
    out={'witness':{'vertices':wp,'spectrum':spectrum(wp)},'central_winding_two_obstruction':{'vertices':ob,'spectrum':spectrum(ob)}}
    out['verdict']='PASS' if out['witness']['spectrum']['winding_values']==[0,1] and 2 in out['central_winding_two_obstruction']['spectrum']['winding_values'] else 'FAIL'
    p=ROOT/ns.json_out;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(jsonable(out),indent=2)+'\n');print(json.dumps(jsonable(out),indent=2));return 0 if out['verdict']=='PASS' else 1
if __name__=='__main__':raise SystemExit(main())
