#!/usr/bin/env python3
r"""
certify_avgdet_neg.py -- RIGOROUS (validated-interval) certificate that the
SEGMENT-AVERAGED Jacobian determinant of the witness is NEGATIVE at the shipped
collision-candidate pair (p,q):   det avg-DF(p,q) < 0,   where

    avg-DF(p,q) = \int_0^1 DF(p + t(q-p)) dt = A * diag(dbar) * B,
    dbar_i = \int_0^1 sigma'(beta_i + alpha_i t) dt
           = (sigma(beta_i + alpha_i) - sigma(beta_i)) / alpha_i     (alpha_i != 0),
           = sigma'(beta_i)                                          (alpha_i = 0),

with beta_i = b_i.p + c_i,  alpha_i = b_i.(q-p),  and sigma the logistic function.
The closed form makes each dbar_i EXACT up to the transcendental sigma, so NO
numerical integration is needed; then det avg-DF = sum_I m_I prod_{i in I} dbar_i
(Cauchy-Binet), m_I = det(A_cols_I) det(B_rows_I).  We enclose each dbar_i with
mpmath.iv and certify that the UPPER endpoint of the det enclosure is < 0.

This upgrades the paper's separation datum (Section 6) from a double-precision
diagnostic to a certified claim.  Combined with det DF(p) > 0 (the global
determinant certificate, at the coincident pair avg-DF(p,p)=DF(p)), the intermediate
value theorem then gives an interior segment fraction at which the averaged Jacobian
is exactly singular.

Run:  python3 code/certify_avgdet_neg.py
Exit 0 iff det avg-DF(p,q) is certified < 0.  No binary float on the decision path
(inputs are exact rationals; sigma is enclosed outward).
"""
from __future__ import annotations
import hashlib, itertools, json, os, sys
from decimal import Decimal
from fractions import Fraction as Fr
from mpmath import iv, mp


def q_to_iv(q: Fr):
    """Rigorous iv enclosure of a Fraction (outward; pad dominates the rounding)."""
    num, den = q.numerator, q.denominator
    approx = mp.mpf(num) / mp.mpf(den)
    ulp = mp.mpf(2) ** (1 - mp.prec)
    pad = abs(approx) * ulp * 8 + mp.mpf(2) ** (4 - mp.prec)
    return iv.mpf([approx - pad, approx + pad])


def sigma_iv(T):
    """logistic sigma(t) = 1/(1 + e^{-t}); denominator >= 1 > 0, reciprocal safe."""
    return 1 / (1 + iv.exp(-T))


def sigmap_iv(T):
    """sigma'(t) = 1/(2 + e^t + e^{-t})."""
    return 1 / (2 + iv.exp(T) + iv.exp(-T))


def _mid(x):
    return (mp.mpf(x.a) + mp.mpf(x.b)) / 2


def main():
    mp.dps = 40
    iv.dps = 40
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, os.pardir, "data", "witness.json")
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    data = json.loads(raw.decode("utf-8"), parse_float=Decimal)
    A = [[Fr(v) for v in row] for row in data["A"]]     # 2x4
    B = [[Fr(v) for v in row] for row in data["B"]]     # 4x2
    c = [Fr(v) for v in data["c"]]                      # 4
    p = [Fr(v) for v in data["witness_pair"]["p"]]      # 2
    q = [Fr(v) for v in data["witness_pair"]["q"]]      # 2
    print(f"witness.json sha256 = {sha}")

    # exact affine data of the segment:  z_i(t) = beta_i + alpha_i t
    dq = [q[0] - p[0], q[1] - p[1]]
    beta = [B[i][0] * p[0] + B[i][1] * p[1] + c[i] for i in range(4)]
    alpha = [B[i][0] * dq[0] + B[i][1] * dq[1] for i in range(4)]

    # dbar_i  (averaged sigma' over the segment; iv)
    dbar = []
    for i in range(4):
        if alpha[i] == 0:
            dbar.append(sigmap_iv(q_to_iv(beta[i])))
        else:
            bi, ai = q_to_iv(beta[i]), q_to_iv(alpha[i])
            dbar.append((sigma_iv(bi + ai) - sigma_iv(bi)) / ai)

    # det avg-DF = sum_I m_I prod_{i in I} dbar_i   (Cauchy-Binet, m_I exact)
    def detAcol(i, j):
        return A[0][i] * A[1][j] - A[0][j] * A[1][i]
    def detBrow(i, j):
        return B[i][0] * B[j][1] - B[i][1] * B[j][0]
    det = iv.mpf(0)
    for I in itertools.combinations(range(4), 2):
        mI = detAcol(*I) * detBrow(*I)
        det = det + q_to_iv(mI) * dbar[I[0]] * dbar[I[1]]

    # avg-DF matrix (2x2) + a diagnostic sigma_min / condition number (from midpoints)
    M = [[iv.mpf(0), iv.mpf(0)], [iv.mpf(0), iv.mpf(0)]]
    for r in range(2):
        for cc in range(2):
            s = iv.mpf(0)
            for i in range(4):
                s = s + q_to_iv(A[r][i]) * dbar[i] * q_to_iv(B[i][cc])
            M[r][cc] = s
    a, b = _mid(M[0][0]), _mid(M[0][1])
    cc_, d = _mid(M[1][0]), _mid(M[1][1])
    fro2 = a*a + b*b + cc_*cc_ + d*d
    dm = a*d - b*cc_
    smax = mp.sqrt((fro2 + mp.sqrt(abs(fro2*fro2 - 4*dm*dm))) / 2)
    smin = abs(dm) / smax

    print(f"avg-DF(p,q) ~ [[{mp.nstr(a,9)}, {mp.nstr(b,9)}],")
    print(f"               [{mp.nstr(cc_,9)}, {mp.nstr(d,9)}]]   (midpoints, diagnostic)")
    print(f"sigma_min ~ {mp.nstr(smin,4)}   kappa_2 ~ {mp.nstr(smax/smin,4)}   (diagnostic)")
    print(f"det avg-DF enclosure = [{mp.nstr(mp.mpf(det.a),8)}, {mp.nstr(mp.mpf(det.b),8)}]")
    print("=" * 60)
    if mp.mpf(det.b) < 0:
        print("VERDICT: PASS -- det avg-DF(p,q) < 0 certified (upper endpoint < 0).")
        return 0
    print("VERDICT: FAIL -- det avg-DF enclosure does not exclude 0 from above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
