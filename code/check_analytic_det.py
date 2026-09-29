#!/usr/bin/env python3
"""Independent check of prop:analytic-det (P4/Q2).
Re-derives alpha,beta,delta by direct 2x2 solve (Fractions), verifies every
algebraic hypothesis exactly, then evaluates the transcendental bound in Arb
(python-flint), a library disjoint from the repo's mpmath.iv path.
Also independently re-proves the proposition's logic numerically by random
sampling of (z2,z3) to confirm the chain of inequalities is not vacuous.
"""

# Assertions below are part of verification, not optional diagnostics.
import sys as _verification_sys
if _verification_sys.flags.optimize:
    raise SystemExit("FAIL: verification requires assertions; omit -O/-OO and unset PYTHONOPTIMIZE.")

import json
from fractions import Fraction as Fr
from pathlib import Path
from report_validation import validate_planar_four
from flint import arb, ctx
ctx.prec = 400

ROOT = Path(__file__).resolve().parents[1]
d = json.loads((ROOT / 'data' / 'witness.json').read_text(encoding='utf-8'), parse_float=Fr, parse_int=Fr)
A = [[Fr(x) for x in r] for r in d['A']]
B = [[Fr(x) for x in r] for r in d['B']]
c = [Fr(x) for x in d['c']]
validate_planar_four(A, B, c)
det2 = lambda u, v: u[0]*v[1] - u[1]*v[0]
Ac = [(A[0][i], A[1][i]) for i in range(4)]
m = {(i,j): det2(Ac[i], Ac[j]) * det2(B[i], B[j]) for i in range(4) for j in range(i+1,4)}

# solve b0 = alpha b2 + beta b3 directly (Gaussian elimination on 2x2)
b0, b2, b3 = B[0], B[2], B[3]
D = det2(b2, b3)
assert D != 0
alpha = det2(b0, b3) / D
beta  = det2(b2, b0) / D
# verify the decomposition EXACTLY
assert alpha*b2[0] + beta*b3[0] == b0[0] and alpha*b2[1] + beta*b3[1] == b0[1], "decomposition fails"
delta = c[0] - alpha*c[2] - beta*c[3]
a, b = abs(alpha), abs(beta)
r = Fr(1) - a - b
mu, lam = Fr(51,1000), Fr(949,1000)

print("alpha ~", float(alpha), " beta ~", float(beta), " delta ~", float(delta))
print("a=|alpha| ~", float(a), " b=|beta| ~", float(b), " r=1-a-b ~", float(r))
H = {
 'r>0': r > 0, 'mu>=a': mu >= a, 'lam>=b': lam >= b, 'mu+lam=1': mu+lam == 1,
 'm02>0': m[(0,2)] > 0, 'm03>0': m[(0,3)] > 0, 'm23<0': m[(2,3)] < 0,
 'm01>=0': m[(0,1)] >= 0, 'm12>=0': m[(1,2)] >= 0, 'm13>=0': m[(1,3)] >= 0,
}
print("exact hypotheses:", H, " ALL:", all(H.values()))
assert all(H.values()), "analytic determinant hypotheses failed"

arbf = lambda fr: arb(fr.numerator)/arb(fr.denominator)
def psi(t):  # 2 log cosh(t/2), evaluated stably
    return 2*(t/2).cosh().log()

t = arbf(abs(delta)/r)
C = arbf(r)*psi(t)
lower = (arbf(m[(0,2)])/arbf(lam))**arbf(lam) * (arbf(m[(0,3)])/arbf(mu))**arbf(mu) * (-C).exp()
neg   = arbf(-m[(2,3)])
margin = lower - neg
print("\nARB dominance lower :", lower.str(25, radius=True))
print("ARB |m23|           :", neg.str(25, radius=True))
print("ARB margin          :", margin.str(25, radius=True))
print("margin > 0 (rigorous):", margin > 0)
print("paper constants hold: lower>0.1290339743302:", lower > arb('0.1290339743302'),
      " |m23|<0.1261291706757:", neg < arb('0.1261291706757'),
      " margin>0.0029048036545:", margin > arb('0.0029048036545'))
assert margin > 0, "analytic determinant margin is not positive"
assert lower > arb('0.1290339743302')
assert neg < arb('0.1261291706757')
assert margin > arb('0.0029048036545')

# sanity: the full chain det/(rho2 rho3) >= m23 + lower at random (z2,z3)
import random
random.seed(7)
sigp = lambda z: 1/(2 + z.exp() + (-z).exp())
worst = None
for _ in range(4000):
    z2 = arb(random.uniform(-60,60)); z3 = arb(random.uniform(-60,60))
    z0 = arbf(alpha)*z2 + arbf(beta)*z3 + arbf(delta)
    # z1 from x: need x; z2,z3 free <-> x = Binv...  z1 = affine(z2,z3):
    # b1 = g2 b2 + g3 b3
    g2 = det2(B[1], b3)/D; g3 = det2(b2, B[1])/D
    z1 = arbf(g2)*z2 + arbf(g3)*z3 + arbf(c[1] - g2*c[2] - g3*c[3])
    rho = {0: sigp(z0), 1: sigp(z1), 2: sigp(z2), 3: sigp(z3)}
    detdf = sum(arbf(m[(i,j)])*rho[i]*rho[j] for (i,j) in m)
    lhs = detdf/(rho[2]*rho[3])
    rhs = arbf(m[(2,3)]) + lower
    ok = lhs >= rhs or float((lhs-rhs).mid()) > -1e-25
    if not ok:
        worst = (float(z2.mid()), float(z3.mid()), float((lhs-rhs).mid()))
print("\nrandom (z2,z3) sanity of det/(r2 r3) >= m23 + LB:", "no violation" if worst is None else worst)
assert worst is None, "sampled analytic inequality violation"
print("VERDICT: INDEPENDENT ANALYTIC DETERMINANT CHECK PASSED")
