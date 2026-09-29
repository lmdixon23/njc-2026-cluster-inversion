#!/usr/bin/env python3
"""Exact verification for the three comparison polygons used in the paper."""
from __future__ import annotations

# Assertions below are part of verification, not optional diagnostics.
import sys as _verification_sys
if _verification_sys.flags.optimize:
    raise SystemExit("FAIL: verification requires assertions; omit -O/-OO and unset PYTHONOPTIMIZE.")

import json
from fractions import Fraction as F
from pathlib import Path
from asymptotic_polygon import add, polygon_simplicity, strict_convexity_by_turns, polygon_signed_double_area, jsonable

def polygon_from_half(g):
    edges=list(g)+[(-x,-y) for x,y in g]
    v=[(F(0),F(0))]
    for e in edges: v.append(add(v[-1],e))
    assert v[-1]==(0,0)
    return edges,v[:-1]

def main():
    # Strictly angularly ordered first-half generators -> convex zonogon boundary.
    convex=[(F(4),F(-3)),(F(4),F(-1)),(F(4),F(1)),(F(4),F(3))]
    simple=[(F(1),F(-3)),(F(1),F(1)),(F(1),F(0)),(F(1),F(3))]
    crossing=[(F(1),F(-3)),(F(7),F(7)),(F(7),F(0)),(F(1),F(3))]
    out={}
    for name,g in [('convex_sign_coherent',convex),('nonconvex_simple',simple),('same_sign_data_self_intersecting',crossing)]:
        e,v=polygon_from_half(g); simp,fail=polygon_simplicity(v); conv,ori=strict_convexity_by_turns(e)
        out[name]={'half_edges':g,'vertices':v,'simple':simp,'strictly_convex':conv,'orientation':ori,'signed_area':polygon_signed_double_area(v)/2,'failures':fail}
    # Positive rescaling preserves every pairwise determinant sign.
    sign_same=True
    for i in range(4):
        for j in range(i+1,4):
            a=simple[i][0]*simple[j][1]-simple[i][1]*simple[j][0]
            b=crossing[i][0]*crossing[j][1]-crossing[i][1]*crossing[j][0]
            sign_same &= (a>0)-(a<0)==(b>0)-(b<0)
    out['sign_pattern_preserved']=sign_same
    passed=(out['convex_sign_coherent']['simple'] and out['convex_sign_coherent']['strictly_convex'] and out['nonconvex_simple']['simple'] and not out['nonconvex_simple']['strictly_convex'] and not out['same_sign_data_self_intersecting']['simple'] and sign_same)
    out['verdict']='PASS' if passed else 'FAIL'
    root=Path(__file__).resolve().parents[1]
    (root/'verification'/'comparison_examples.json').write_text(json.dumps(jsonable(out),indent=2)+'\n', newline="\n")
    print(json.dumps(jsonable(out),indent=2)); raise SystemExit(0 if passed else 1)
if __name__=='__main__': main()
