#!/usr/bin/env python3
"""From-scratch exact-rational asymptotic polygon toolkit (independent of
asymptotic_polygon.py). Implements:
  * sector construction of the polygon from (A,B) per lem:polygon/prop:normal-form
  * exact simplicity via segment predicates (all pairs)
  * exact winding-number function and full winding SPECTRUM via arrangement
    face sampling (fragment midpoints offset by certified-small epsilon)
All arithmetic in fractions.Fraction.
"""
from fractions import Fraction as Fr
from itertools import combinations

def det2(u, v): return u[0]*v[1] - u[1]*v[0]
def sub(u, v): return (u[0]-v[0], u[1]-v[1])
def add(u, v): return (u[0]+v[0], u[1]+v[1])

# ---------- polygon from (A,B) ----------
def half_plane_angle_key(v):
    """Exact angular sort key for nonzero rational vectors over [0,2pi)."""
    x, y = v
    if y > 0 or (y == 0 and x > 0):
        upper = 0  # angle in [0,pi)
    else:
        upper = 1
    return (upper, )  # combined with cross-product comparator below

def angular_sort(dirs):
    """Sort rational direction vectors by angle in [0,2pi), exactly."""
    def half(v):
        x, y = v
        return 0 if (y > 0 or (y == 0 and x > 0)) else 1
    import functools
    def cmp(u, v):
        hu, hv = half(u[0]), half(v[0])
        if hu != hv: return -1 if hu < hv else 1
        d = det2(u[0], v[0])
        if d > 0: return -1
        if d < 0: return 1
        return 0
    return sorted(dirs, key=functools.cmp_to_key(cmp))

def build_polygon_indep(A, B):
    """Vertices of the asymptotic polygon, ccw by angular traversal of u."""
    N = len(B)
    acols = [(A[0][i], A[1][i]) for i in range(N)]
    walls = []  # (direction w with b_i . w = 0, row i)
    for i, b in enumerate(B):
        w = (-b[1], b[0])
        walls.append((w, i)); walls.append(((b[1], -b[0]), i))
    walls = angular_sort(walls)
    K = len(walls)
    verts = []
    crossings = []
    for k in range(K):
        d1 = walls[k][0]; d2 = walls[(k+1) % K][0]
        u = add(d1, d2)  # interior direction of sector (gap < pi by antipodality)
        assert det2(d1, u) > 0 and det2(u, d2) > 0, "u not interior"
        s = [1 if (B[i][0]*u[0] + B[i][1]*u[1]) > 0 else 0 for i in range(N)]
        assert all(B[i][0]*u[0] + B[i][1]*u[1] != 0 for i in range(N))
        v = (sum(acols[i][0] for i in range(N) if s[i]),
             sum(acols[i][1] for i in range(N) if s[i]))
        verts.append(v)
        crossings.append(walls[(k+1) % K][1])
    return verts, crossings, acols

# ---------- exact segment predicates ----------
def on_seg(x, p, q):
    return det2(sub(q,p), sub(x,p)) == 0 and \
           (x[0]-p[0])*(x[0]-q[0]) <= 0 and (x[1]-p[1])*(x[1]-q[1]) <= 0

def segs_meet(p1, q1, p2, q2):
    o1 = det2(sub(q1,p1), sub(p2,p1)); o2 = det2(sub(q1,p1), sub(q2,p1))
    o3 = det2(sub(q2,p2), sub(p1,p2)); o4 = det2(sub(q2,p2), sub(q1,p2))
    if ((o1>0) != (o2>0)) and ((o3>0) != (o4>0)) and o1 != 0 and o2 != 0 and o3 != 0 and o4 != 0:
        return True
    return on_seg(p2,p1,q1) or on_seg(q2,p1,q1) or on_seg(p1,p2,q2) or on_seg(q1,p2,q2)

def is_simple(verts):
    n = len(verts)
    if len(set(verts)) != n: return False, 'repeated vertex'
    E = [(verts[i], verts[(i+1)%n]) for i in range(n)]
    for (p,q) in E:
        if p == q: return False, 'zero edge'
    for i in range(n):
        j = (i+1) % n
        # adjacent: shared endpoint only
        if on_seg(E[i][0], *E[j]) or on_seg(E[j][1], *E[i]):
            return False, f'adjacent overlap {i},{j}'
    for i, j in combinations(range(n), 2):
        if j == (i+1)%n or i == (j+1)%n: continue
        if segs_meet(*E[i], *E[j]):
            return False, f'nonadjacent meet {i},{j}'
    return True, 'ok'

# ---------- winding ----------
def winding(verts, y):
    """Exact winding number of the closed polyline around y (y not on curve)."""
    n = len(verts); w = 0
    for i in range(n):
        p, q = verts[i], verts[(i+1)%n]
        assert not on_seg(y, p, q), 'query on curve'
        if (p[1] <= y[1]) and (q[1] > y[1]):
            if det2(sub(q,p), sub(y,p)) > 0: w += 1
        elif (p[1] > y[1]) and (q[1] <= y[1]):
            if det2(sub(q,p), sub(y,p)) < 0: w -= 1
    return w

def seg_intersection_point(p1, q1, p2, q2):
    d1, d2 = sub(q1,p1), sub(q2,p2)
    den = det2(d1, d2)
    if den == 0: return None
    t = det2(sub(p2,p1), d2) / den
    s = det2(sub(p2,p1), d1) / den
    if 0 <= t <= 1 and 0 <= s <= 1:
        return (p1[0]+t*d1[0], p1[1]+t*d1[1])
    return None

def winding_spectrum(verts):
    """Exact spectrum over all complementary faces: cut edges at all crossing
    points, sample each fragment midpoint offset both ways by an epsilon
    certified smaller than the distance to every non-incident fragment."""
    n = len(verts)
    E = [(verts[i], verts[(i+1)%n]) for i in range(n)]
    cuts = [set([Fr(0), Fr(1)]) for _ in range(n)]
    for i, j in combinations(range(n), 2):
        p1,q1 = E[i]; p2,q2 = E[j]
        d1, d2 = sub(q1,p1), sub(q2,p2)
        den = det2(d1,d2)
        if den == 0: continue
        t = det2(sub(p2,p1), d2)/den; s = det2(sub(p2,p1), d1)/den
        if 0 <= t <= 1 and 0 <= s <= 1:
            cuts[i].add(t); cuts[j].add(s)
    samples = []
    frags = []
    for i in range(n):
        p, q = E[i]; d = sub(q,p)
        ts = sorted(cuts[i])
        for a, b in zip(ts, ts[1:]):
            tm = (a+b)/2
            mid = (p[0]+tm*d[0], p[1]+tm*d[1])
            nrm = (-d[1], d[0])
            frags.append((p, q, a, b))
            samples.append((mid, nrm))
    def dist2_point_seg(x, p, q):
        d = sub(q,p); L2 = d[0]*d[0]+d[1]*d[1]
        t = ( (x[0]-p[0])*d[0] + (x[1]-p[1])*d[1] ) / L2
        t = max(Fr(0), min(Fr(1), t))
        c = (p[0]+t*d[0], p[1]+t*d[1])
        dx, dy = x[0]-c[0], x[1]-c[1]
        return dx*dx+dy*dy
    spec = set()
    for mid, nrm in samples:
        L2 = nrm[0]*nrm[0]+nrm[1]*nrm[1]
        dmin2 = None
        for (p,q) in E:
            if on_seg(mid, p, q):   # incident (mid lies on this edge)
                continue
            dd = dist2_point_seg(mid, p, q)
            dmin2 = dd if dmin2 is None else min(dmin2, dd)
        # eps so that |eps * nrm| < sqrt(dmin2)/2   -> eps^2 L2 < dmin2/4
        if dmin2 is None or dmin2 == 0:
            eps = Fr(1, 10**12)
        else:
            eps = Fr(1, 4)
            while eps*eps*L2 >= dmin2/4:
                eps /= 2
        for sgn in (1, -1):
            y = (mid[0]+sgn*eps*nrm[0], mid[1]+sgn*eps*nrm[1])
            if any(on_seg(y, p, q) for (p,q) in E):
                continue
            spec.add(winding(verts, y))
    # far outside always 0
    spec.add(0)
    return sorted(spec)

def poly_from_half_edges(half):
    edges = list(half) + [(-e[0], -e[1]) for e in half]
    v = [(Fr(0), Fr(0))]
    for e in edges:
        v.append(add(v[-1], e))
    assert v[-1] == v[0]
    return v[:-1]

if __name__ == '__main__':
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    d = json.loads((root / 'data' / 'witness.json').read_text(encoding='utf-8'), parse_float=Fr, parse_int=Fr)
    A = [[Fr(x) for x in r] for r in d['A']]; B = [[Fr(x) for x in r] for r in d['B']]
    verts, crossings, acols = build_polygon_indep(A, B)
    print('independent crossing row order (per sector step):', crossings)
    simple, why = is_simple(verts)
    print('independent simplicity:', simple, why)
    # signed double area
    n = len(verts)
    dbl = sum(det2(verts[i], verts[(i+1)%n]) for i in range(n))
    print('signed area exact =', dbl/2, '~', float(dbl/2))
    spec = winding_spectrum(verts)
    print('witness winding spectrum =', spec)
    # interior point check q=(-25,4)? that's the FUNCTIONAL, not a point. Use centroid-ish.
    # sign-insufficiency polygons (prop:sign-insufficient)
    P1 = poly_from_half_edges([(Fr(1),Fr(-3)),(Fr(1),Fr(1)),(Fr(1),Fr(0)),(Fr(1),Fr(3))])
    P2 = poly_from_half_edges([(Fr(1),Fr(-3)),(Fr(7),Fr(7)),(Fr(7),Fr(0)),(Fr(1),Fr(3))])
    print('\nP1 vertices:', [(str(a),str(b)) for a,b in P1])
    print('P1 simple:', is_simple(P1))
    print('P2 vertices:', [(str(a),str(b)) for a,b in P2])
    print('P2 simple:', is_simple(P2))
    # exact crossing points of P2
    E = [(P2[i], P2[(i+1)%8]) for i in range(8)]
    xs = []
    for i, j in combinations(range(8), 2):
        if j == (i+1)%8 or i == (j+1)%8: continue
        pt = seg_intersection_point(*E[i], *E[j])
        if pt is not None:
            xs.append(((i,j), (str(pt[0]), str(pt[1]))))
    print('P2 nonadjacent crossings:', xs)
    print('P2 winding spectrum:', winding_spectrum(P2))
    # matching pairwise determinant signs between P1 and P2 edge lists
    H1 = [(Fr(1),Fr(-3)),(Fr(1),Fr(1)),(Fr(1),Fr(0)),(Fr(1),Fr(3))]
    H2 = [(Fr(1),Fr(-3)),(Fr(7),Fr(7)),(Fr(7),Fr(0)),(Fr(1),Fr(3))]
    sgn = lambda x: (x>0)-(x<0)
    same = all(sgn(det2(H1[i],H1[j]))==sgn(det2(H2[i],H2[j])) for i,j in combinations(range(4),2))
    print('pairwise det signs identical:', same)
    if not simple or spec != [0, 1] or not same:
        raise SystemExit('VERDICT: INDEPENDENT POLYGON CHECK FAILED')
    print('VERDICT: INDEPENDENT POLYGON CHECK PASSED')
