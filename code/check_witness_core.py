#!/usr/bin/env python3
"""Independent re-derivation (Q6 core + thm:separation numbers).
Library: fractions.Fraction (exact) + python-flint arb (balls), NOT the repo's
mpmath.iv / decimal_interval code. Formulas re-derived from the manuscript:
  det DF(x) = sum_{i<j} m_ij s'(z_i) s'(z_j),  m_ij = det(A_{:,ij}) det(B_{ij,:})
  avg DF over segment: entrywise integral -> diag entries are divided
  differences (sigma(z1)-sigma(z0))/(z1-z0) since z_i is affine in t.
"""
import json, hashlib
from fractions import Fraction as Fr
from pathlib import Path
from flint import arb, ctx

ctx.prec = 200  # bits

ROOT = Path(__file__).resolve().parents[1]
raw = (ROOT / 'data' / 'witness.json').read_bytes()
witness_sha = hashlib.sha256(raw).hexdigest()
print('file sha256 =', witness_sha)
d = json.loads(raw, parse_float=Fr, parse_int=Fr)

A = [[Fr(x) for x in row] for row in d['A']]   # Fraction(float) is exact
B = [[Fr(x) for x in row] for row in d['B']]
c = [Fr(x) for x in d['c']]
p = [Fr(x) for x in d['witness_pair']['p']]
q = [Fr(x) for x in d['witness_pair']['q']]

def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]

# --- exact mixed minors ---
Acols = [[A[0][j], A[1][j]] for j in range(4)]
Brows = B
m = {}
for i in range(4):
    for j in range(i+1, 4):
        m[(i,j)] = det2(Acols[i], Acols[j]) * det2(Brows[i], Brows[j])

claimed = d['claimed']['minors']
print('\n--- exact mixed minors (Fraction) vs claimed floats ---')
ok_all = True
for (i,j), val in m.items():
    key = f'{i}{j}'
    cl = claimed[key]
    ok = abs(float(val) - float(cl)) <= 1e-12
    ok_all &= ok
    print(f'm_{i}{j} exact~{float(val):+.17g}  claimed {float(cl):+.17g}  display_match={ok}')
negative_display_match = abs(float(m[(2,3)]) - float(d['claimed']['negative_minor_value'])) <= 1e-12
print('claimed negative_minor_value display match:', negative_display_match)
signs = sorted([(i,j) for (i,j),v in m.items() if v > 0]), sorted([(i,j) for (i,j),v in m.items() if v < 0])
print('positive pairs:', signs[0], ' negative pairs:', signs[1])

# Bind the independent exact signs to the regenerated point-tail report.  The
# current report deliberately records signs rather than duplicating interval
# strings from the generating implementation.
ptr = json.loads((ROOT / 'results' / 'point_tail_recomputed.json').read_text(encoding='utf-8'))
expected_positive = [list(pair) for pair in signs[0]]
expected_negative = [list(pair) for pair in signs[1]]
point_report_bound = (
    ptr.get('witness_sha256') == witness_sha
    and ptr.get('minor_signs') == {
        'positive': expected_positive,
        'negative': expected_negative,
        'ambiguous': [],
    }
)
print('point-tail witness and exact minor signs agree:', point_report_bound)

# --- z values at p, q (exact) ---
def z_at(x):
    return [B[i][0]*x[0] + B[i][1]*x[1] + c[i] for i in range(4)]
zp, zq = z_at(p), z_at(q)

def arbf(fr):  # exact Fraction -> arb ball
    return arb(fr.numerator) / arb(fr.denominator)

def sig(t):   # logistic in arb
    return 1/(1+(-t).exp())
def sigp(t):  # sigma' = 1/(2+e^t+e^-t)
    return 1/(2 + t.exp() + (-t).exp())

# --- averaged Jacobian determinant, closed form ---
avg = []
for i in range(4):
    z0, z1 = arbf(zp[i]), arbf(zq[i])
    dz = arbf(zq[i]-zp[i])
    if zp[i] == zq[i]:
        avg.append(sigp(z0))
    else:
        avg.append((sig(z1)-sig(z0))/dz)

detavg = arb(0)
for (i,j), mij in m.items():
    detavg += arbf(mij)*avg[i]*avg[j]
print('\n--- averaged Jacobian ---')
print('det avg DF ball =', detavg.str(30, radius=True))
print('claimed         =', float(d['claimed']['seg_avg_det']))
print('certified negative:', detavg < 0)
claimed_avg_close = abs(float(detavg.mid()) - float(d['claimed']['seg_avg_det'])) <= 1e-12
print('claimed display value agrees to 1e-12:', claimed_avg_close)

# also the averaged matrix entries + sigma_min for sec:reconcile numbers
Mbar = [[arb(0),arb(0)],[arb(0),arb(0)]]
for r in range(2):
    for s_ in range(2):
        acc = arb(0)
        for i in range(4):
            acc += arbf(A[r][i])*avg[i]*arbf(B[i][s_])
        Mbar[r][s_] = acc
print('avgDF matrix ~ [[%s, %s],[%s, %s]]' % tuple(Mbar[r][s_].str(10) for r in range(2) for s_ in range(2)))
# singular values via eigen of M^T M (2x2 closed form)
a_,b_,c_,d_ = Mbar[0][0],Mbar[0][1],Mbar[1][0],Mbar[1][1]
t1 = a_*a_+b_*b_+c_*c_+d_*d_
dt = a_*d_-b_*c_
disc = (t1*t1 - 4*dt*dt).sqrt()
smax = ((t1+disc)/2).sqrt(); smin = ((t1-disc)/2).sqrt()
print('sigma_min ~', smin.str(8), ' kappa2 ~', (smax/smin).str(8))

# --- min det DF on the segment: scan + refine (mpmath-free, arb) ---
def detDF_at_t(t):  # t arb in [0,1]
    zs = [arbf(zp[i]) + t*arbf(zq[i]-zp[i]) for i in range(4)]
    rho = [sigp(zz) for zz in zs]
    s = arb(0)
    for (i,j), mij in m.items():
        s += arbf(mij)*rho[i]*rho[j]
    return s

N = 4000
vals = []
for k in range(N+1):
    t = arb(k)/N
    vals.append((float(detDF_at_t(t).mid()), k))
vals.sort()
print('\n--- min det DF on segment (scan of %d pts) ---' % (N+1))
kbest = vals[0][1]
# golden-ish refine around kbest
lo = max(0, kbest-1)/N; hi = min(N, kbest+1)/N
for _ in range(60):
    t1 = lo + (hi-lo)/3; t2 = hi - (hi-lo)/3
    if float(detDF_at_t(arb(t1)).mid()) < float(detDF_at_t(arb(t2)).mid()):
        hi = t2
    else:
        lo = t1
tm = (lo+hi)/2
vmin = detDF_at_t(arb(tm))
print('argmin t ~ %.10f   min det DF ~ %s' % (tm, vmin.str(20)))
print('claimed sampled_min_detDF_on_segment (201-point sample, not a certified lower bound) =', float(d['claimed']['min_detDF_on_segment']))
print('scan-min positive:', vals[0][0] > 0)
print('\nall minors double-match:', ok_all)
if not ok_all or not negative_display_match:
    raise SystemExit('VERDICT: INDEPENDENT WITNESS CHECK FAILED: mixed-minor display values')
if not point_report_bound:
    raise SystemExit('VERDICT: INDEPENDENT WITNESS CHECK FAILED: point-tail report binding')
if not detavg < arb('-1.7199e-5'):
    raise SystemExit('VERDICT: INDEPENDENT WITNESS CHECK FAILED: theorem bound -1.7199e-5')
print('theorem bound det average < -1.7199e-5: certified')
if not claimed_avg_close:
    raise SystemExit('VERDICT: INDEPENDENT WITNESS CHECK FAILED: averaged-determinant display value')
if vals[0][0] <= 0:
    raise SystemExit('VERDICT: INDEPENDENT WITNESS CHECK FAILED: segment scan')
print('VERDICT: INDEPENDENT WITNESS CHECK PASSED')
