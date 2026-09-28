#!/usr/bin/env python3
r"""
certify_det_exact.py -- RIGOROUS (validated-interval) certificate that
det DF(x) > 0 for all x in R^2, for F(x)=A*sigma(Bx+c), read DIRECTLY from the
frozen data/witness.json.

WHY: an earlier prototype certificate established det DF>0 using IEEE-double
("rigorous modulo rounding"), read an untracked pickle (not the sha'd
witness.json), computed no sha, and did not fail-exit on unresolved boxes.
This module removes those gaps for the injectivity corollary:
  * inputs are the exact rationals of witness.json (no pickle, no binary-float
    contamination of A,B,c or the minor products);
  * sigma'(t) = 1/(2 + 2 cosh t) is enclosed with mpmath.iv (outward-rounded);
  * det DF over a box is enclosed by an interval; a box is certified only if the
    LOWER end of that interval is > 0;
  * the tail lemma (domination outside radius R0) is reproduced with a rigorous
    outward psi<0 test + Lipschitz slack (same analytic domination mechanism as the earlier prototype);
  * ANY unresolved box (below HMIN) or budget exhaustion => print FAIL, exit 1.

Structure of det DF (Cauchy-Binet):
  det DF(x) = sum_{I} m_I * prod_{i in I} sigma'(t_i),  t_i = b_i . x + c_i,
  m_I = det(A_columns_I) * det(B_rows_I).  Exactly one m_I is negative here: I=(2,3).

Run:
  python3 code/certify_det_exact.py --selftest        # fast primitive checks
  python3 code/certify_det_exact.py                    # selftest + full certificate
  python3 code/certify_det_exact.py --dps 40 --hmin 0.008 --budget 40000000
Exit 0 iff selftest passes AND the tail lemma holds AND every interior box is
certified (det lower bound > 0). Exit 1 on any failure. Interval math is SOUND at
any --dps (lower dps only widens enclosures -> more subdivision, never a false PASS).
Runtime note: mpmath is not vectorized; the interior B&B may take minutes. Progress
is printed. If it is too slow, raise --hmin first (coarser floor) then lower --dps.
"""
from __future__ import annotations
import argparse, json, os, sys, time, itertools
from decimal import Decimal
from fractions import Fraction as Fr
import mpmath
from mpmath import iv, mp
from report_validation import validate_planar_four, require


# --------------------------------------------------------------------------
# exact input
# --------------------------------------------------------------------------
def load_witness(path):
    raw = open(path, "rb").read()
    import hashlib
    sha = hashlib.sha256(raw).hexdigest()
    data = json.loads(raw.decode("utf-8"), parse_float=Decimal)
    A = [[Fr(v) for v in row] for row in data["A"]]   # 2x4
    B = [[Fr(v) for v in row] for row in data["B"]]   # 4x2
    c = [Fr(v) for v in data["c"]]                    # 4
    validate_planar_four(A, B, c)
    return A, B, c, sha


def minor_products(A, B):
    """m_I = det(A_cols_I)*det(B_rows_I), exact Fraction, for all 2-subsets I."""
    validate_planar_four(A, B)
    m = {}
    for I in itertools.combinations(range(4), 2):
        i, j = I
        dA = A[0][i] * A[1][j] - A[0][j] * A[1][i]
        dB = B[i][0] * B[j][1] - B[i][1] * B[j][0]
        m[I] = dA * dB
    return m


# --------------------------------------------------------------------------
# rigorous rational -> interval (outward, independent of mpmath int rounding)
# --------------------------------------------------------------------------
def q_to_iv(q: Fr):
    """Rigorous iv enclosure of a Python Fraction at the current iv precision.
    approx = mpf(num)/mpf(den) incurs <= ~1.5 ulp total rounding; pad by 8 ulp
    (relative) plus a tiny absolute floor => the true value is enclosed."""
    num, den = q.numerator, q.denominator
    approx = mp.mpf(num) / mp.mpf(den)
    ulp = mp.mpf(2) ** (1 - mp.prec)                 # ~2 * 2^-prec (relative)
    pad = abs(approx) * ulp * 8 + mp.mpf(2) ** (4 - mp.prec)
    return iv.mpf([approx - pad, approx + pad])


def sigma_prime_iv(T):
    """sigma'(t) = 1/(2 + 2 cosh t) = 1/(2 + e^t + e^{-t}), enclosed for t in T.
    mpmath.iv has no cosh, so use exp (which it does have). Denominator
    2 + e^t + e^{-t} >= 4 > 0 (AM-GM), so the reciprocal is safe. SOUND: for the
    exact range check, 2 + exp(T) + exp(-T) contains [2+2cosh(t) : t in T] (the
    interval sum only over-widens via the t-dependency, never under-encloses)."""
    return 1 / (2 + iv.exp(T) + iv.exp(-T))


def iv_abs(X):
    """Rigorous |X| for an iv interval (do not rely on __abs__ being defined)."""
    if X.a >= 0:
        return X
    if X.b <= 0:
        return -X
    return iv.mpf([0, max(-X.a, X.b)])


def _mpf_q(q):
    """mpf from a Fraction via exact numerator/denominator (mp.mpf rejects Fraction)."""
    return mp.mpf(q.numerator) / mp.mpf(q.denominator)


# --------------------------------------------------------------------------
# interior: det DF lower bound over an exact-corner box
# --------------------------------------------------------------------------
def t_interval(bi_x: Fr, bi_y: Fr, ci: Fr, x1, y1, x2, y2):
    """Exact [min,max] of b_i . x + c_i over the rectangle (linear => at corners),
    returned as an iv interval."""
    xlo = x1 if bi_x >= 0 else x2
    xhi = x2 if bi_x >= 0 else x1
    ylo = y1 if bi_y >= 0 else y2
    yhi = y2 if bi_y >= 0 else y1
    tmin = bi_x * xlo + bi_y * ylo + ci
    tmax = bi_x * xhi + bi_y * yhi + ci
    lo = q_to_iv(tmin)
    hi = q_to_iv(tmax)
    return iv.mpf([lo.a, hi.b])


def det_lower_positive(box, B, c, m_iv, pos, neg):
    """Return (certified: bool, lower: mpf) for det DF over `box` = (x1,y1,x2,y2)."""
    x1, y1, x2, y2 = box
    sp = [sigma_prime_iv(t_interval(B[i][0], B[i][1], c[i], x1, y1, x2, y2))
          for i in range(4)]
    acc = iv.mpf(0)
    for J in pos:
        acc = acc + m_iv[J] * sp[J[0]] * sp[J[1]]
    In = neg[0]
    acc = acc - m_iv[In] * sp[In[0]] * sp[In[1]]      # m_iv[In] is |m_neg|>0 interval
    return bool(acc.a > 0), mp.mpf(acc.a)             # coerce endpoint to bare mpf


# --------------------------------------------------------------------------
# tail lemma: det DF>0 for ||x||_2 > R0  (faithful rigorous port of certify.py)
# --------------------------------------------------------------------------
def _bnorm_upper(B, i):
    """Rigorous UPPER bound on ||b_i|| (mp.sqrt rounds to nearest; pad up)."""
    return mp.sqrt(_mpf_q(B[i][0])**2 + _mpf_q(B[i][1])**2) * (1 + mp.mpf(2)**(-20))


def _tail_setup(B, c, m, pos, neg):
    """Shared tail constants: rho(u)=|b_i.u| (iv), K_J (iv, outward), and the
    per-pair Lipschitz constant L_J of g_J(u)=rho_neg(u)-rho_J(u)."""
    In = neg[0]
    mneg = -m[In]                                    # >0
    bnf = [(Fr(B[i][0]), Fr(B[i][1])) for i in range(4)]
    def rho(iu, cti, sti):
        return iv_abs(q_to_iv(bnf[iu][0]) * cti + q_to_iv(bnf[iu][1]) * sti)
    KJ = {}
    for J in pos:
        KJ[J] = (iv.log(16 * q_to_iv(mneg) / q_to_iv(m[J]))
                 + q_to_iv(abs(c[J[0]])) + q_to_iv(abs(c[J[1]]))
                 + q_to_iv(abs(c[In[0]])) + q_to_iv(abs(c[In[1]])))
    LJ = {J: (_bnorm_upper(B, In[0]) + _bnorm_upper(B, In[1])
              + _bnorm_upper(B, J[0]) + _bnorm_upper(B, J[1])) for J in pos}
    return In, rho, KJ, LJ


def tail_certify_fixed_R(R, B, c, m, pos, neg, Mgrid=200000, verbose=False):
    """RIGOROUS: certify det DF(x) > 0 for all ||x||_2 > R (R a Python int).

    Exact sufficient condition (from  sum|t_neg| - sum|t_J| >= r*g_J - sum|c|  and
    sigma'(t) in [(1/4)e^-|t|, e^-|t|]): for every direction u some POSITIVE pair J
    dominates the single negative term at radius R,
        R * g_J(u) > K_J     with    g_J(u) = rho_neg(u) - rho_J(u) > 0,
    rho_i(u)=|b_i.u|,  K_J = ln(16|m_neg|/m_J) + |c_2|+|c_3|+|c_j0|+|c_j1|.
    Between grid points g_J is controlled by its Lipschitz constant
    L_J = ||b_2||+||b_3||+||b_j0||+||b_j1|| (each |b_i.u| is ||b_i||-Lipschitz in the
    angle): over an arc of half-width dtheta/2, g_J drops by at most L_J*dtheta/2,
    so g_arc_lo = g_J(theta_k) - L_J*dtheta/2 is a rigorous lower bound over the arc.
    A direction certifies iff some J has g_arc_lo > 0 AND R*g_arc_lo - K_J^upper > 0.
    All bounds are outward, so a True verdict is rigorous; larger R only helps
    (monotone). The dominating pair may SWITCH between directions.

    This REPLACES an earlier flawed step that took sup_theta min_J K_J/g_J and then
    inflated it by min_theta max_J g_J -- combining the ratio-minimizing pair with a
    gap that need not belong to it (separate internal reconstruction, grid index 37576).
    """
    require(type(R) is int and R > 0 and type(Mgrid) is int and Mgrid > 0, 'R and grid must be positive integers')
    In, rho, KJ, LJ = _tail_setup(B, c, m, pos, neg)
    half_arc = mp.pi / Mgrid * (1 + mp.mpf(2)**(-30))    # scalar upper bound on dtheta/2
    dth = 2 * iv.pi / Mgrid
    worst = mp.mpf('+inf')                               # min over grid of best margin
    for k in range(Mgrid):
        th = dth * k
        cti, sti = iv.cos(th), iv.sin(th)
        rhon = rho(In[0], cti, sti) + rho(In[1], cti, sti)
        best = mp.mpf('-inf')
        for J in pos:
            rhoJ = rho(J[0], cti, sti) + rho(J[1], cti, sti)
            # Keep g_J and the margin as INTERVALS; accept only on their lower
            # endpoints. The only scalar is `drop`, an explicit conservative UPPER
            # bound on the arc's Lipschitz decrease of g_J, subtracted as a point.
            drop = LJ[J] * half_arc
            g_iv = (rhon - rhoJ) - iv.mpf(drop)                     # iv; g_iv.a = arc lower bound of g_J
            if not (mp.mpf(g_iv.a) > 0):                            # g_J>0 across the arc (needed for r>=R)
                continue
            margin = mp.mpf((R * g_iv - KJ[J]).a)                   # iv lower bound of R*g_J - K_J
            if margin > best:
                best = margin
        if best < worst:
            worst = best
        if not (best > 0):
            if verbose:
                print(f"    [R={R}] direction k={k} NOT dominated "
                      f"(best margin {mp.nstr(best,4)})")
            return False, worst
    return True, worst


# --------------------------------------------------------------------------
# selftest -- validate every primitive fast, before any long run
# --------------------------------------------------------------------------
def selftest(A, B, c, m, pos, neg, m_iv):
    ok = True
    def chk(name, cond):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        ok = ok and cond

    # (a) sigma' enclosures
    s0 = sigma_prime_iv(iv.mpf([0, 0]))
    chk("sigma'(0) encloses 1/4", s0.a <= mp.mpf(1)/4 <= s0.b)
    for t in ('0.7', '-2.3', '5.0'):
        T = iv.mpf([t, t]); S = sigma_prime_iv(T)
        ref = 1/(2 + 2*mp.cosh(mp.mpf(t)))
        chk(f"sigma'({t}) encloses ref", S.a <= ref <= S.b)
    big = sigma_prime_iv(iv.mpf([20, 40]))
    chk("sigma' tiny for large t", big.b < mp.mpf('1e-6'))

    # (b) exact minor products: exactly one negative, at (2,3)
    negs = [I for I in m if m[I] < 0]
    chk("exactly one negative minor product", len(negs) == 1)
    chk("negative minor product is (2,3)", negs == [(2, 3)])
    # cross-check magnitudes vs polygon checker / paper
    approx = {(0,1):0.0333917,(0,2):0.3887019,(0,3):0.0008636,
              (1,2):3.4206863,(1,3):0.0086087,(2,3):-0.1261292}
    good = all(abs(float(m[I]) - approx[I]) < 1e-6 for I in approx)
    chk("minor products match paper values", good)

    # (c) det lower bound contains the true det at a sample point (degenerate box)
    x0 = (Fr(0), Fr(0), Fr(0), Fr(0))
    cert0, lo0 = det_lower_positive(x0, B, c, m_iv, pos, neg)
    sp_ref = [1/(2 + 2*mp.cosh(mp.mpf(c[i].numerator)/c[i].denominator)) for i in range(4)]
    det0 = mp.mpf(0)
    for J in pos:
        det0 += (mp.mpf(m[J].numerator)/m[J].denominator) * sp_ref[J[0]] * sp_ref[J[1]]
    In = neg[0]
    det0 -= (mp.mpf((-m[In]).numerator)/(-m[In]).denominator) * sp_ref[In[0]] * sp_ref[In[1]]
    chk("det interval lower <= true det at origin", lo0 <= det0)
    chk("origin certifies (det>0 there)", cert0)

    # (d) vacuity guard: with an inflated fake |m_neg| the gate must REFUSE somewhere
    fake = dict(m_iv); fake[neg[0]] = q_to_iv(Fr(1000))
    cert_fake, _ = det_lower_positive((Fr(-1),Fr(-1),Fr(1),Fr(1)), B, c, fake, pos, neg)
    chk("gate is not vacuous (fake huge negative term refuses)", not cert_fake)

    # (e) tail certifier rejects an obviously-too-small radius (fast: early exit)
    okbad, _ = tail_certify_fixed_R(1, B, c, m, pos, neg, Mgrid=4096)
    chk("tail certifier rejects R=1 (too small)", not okbad)
    return ok


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--witness", default=None)
    ap.add_argument("--dps", type=int, default=40)
    ap.add_argument("--hmin", type=float, default=0.008)
    ap.add_argument("--budget", type=int, default=40_000_000)
    ap.add_argument("--grid", type=int, default=200000)
    ap.add_argument("--selftest", action="store_true", help="run only the selftest")
    args = ap.parse_args()
    mp.dps = max(args.dps, 30)
    iv.dps = max(args.dps, 30)

    path = args.witness or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, "data", "witness.json")
    A, B, c, sha = load_witness(path)
    m = minor_products(A, B)
    pos = [I for I in m if m[I] > 0]
    neg = [I for I in m if m[I] < 0]
    m_iv = {I: q_to_iv(m[I] if m[I] > 0 else -m[I]) for I in m}  # store |m_I| as +iv

    print(f"witness.json  path = {os.path.normpath(path)}")
    print(f"witness.json  sha256 = {sha}")
    print(f"dps={mp.dps}  hmin={args.hmin}  budget={args.budget}  grid={args.grid}")
    print("-" * 72)
    print("SELFTEST:")
    if len(neg) != 1:
        print("  [FAIL] expected exactly one negative minor product"); sys.exit(1)
    if not selftest(A, B, c, m, pos, neg, m_iv):
        print("SELFTEST FAILED"); sys.exit(1)
    print("SELFTEST PASSED")
    if args.selftest:
        sys.exit(0)

    print("-" * 72)
    print("TAIL (det DF>0 outside radius R0) -- fixed-R certification:")
    R0 = None
    for Rtry in (582, 583, 585, 600):
        okR, worst = tail_certify_fixed_R(Rtry, B, c, m, pos, neg, Mgrid=args.grid)
        print(f"  R={Rtry}: {'CERTIFIED' if okR else 'not certified'} "
              f"(worst-direction margin {mp.nstr(worst, 4)})")
        if okR:
            R0 = Rtry
            break
    if R0 is None:
        print("  [FAIL] tail NOT certified at tried radii; raise --grid"); sys.exit(1)
    print(f"  [PASS] tail certified at R0 = {R0}")
    R0f = Fr(R0)

    print("-" * 72)
    print(f"INTERIOR branch-and-bound on [-{int(R0)},{int(R0)}]^2 (mpmath.iv):")
    hmin = Fr(args.hmin).limit_denominator(10**9)
    root = (-R0f, -R0f, R0f, R0f)
    stack = [root]
    processed = certified = 0
    fails = []
    t0 = time.time()
    while stack and processed < args.budget:
        box = stack.pop()
        processed += 1
        cert, lo = det_lower_positive(box, B, c, m_iv, pos, neg)
        if cert:
            certified += 1
        else:
            x1, y1, x2, y2 = box
            w, h = x2 - x1, y2 - y1
            if w < hmin and h < hmin:
                fails.append((box, lo))
                if len(fails) <= 3:
                    print(f"  [FAIL] uncertified box < hmin at {[float(v) for v in box]} "
                          f"det_lower={mp.nstr(lo,3)}")
            elif w >= h:
                mx = (x1 + x2) / 2
                stack.append((x1, y1, mx, y2)); stack.append((mx, y1, x2, y2))
            else:
                my = (y1 + y2) / 2
                stack.append((x1, y1, x2, my)); stack.append((x1, my, x2, y2))
        if processed % 500000 == 0:
            print(f"    ... processed={processed} certified={certified} "
                  f"stack={len(stack)} t={time.time()-t0:.0f}s")

    budget_exhausted = bool(stack)
    print("-" * 72)
    print(f"processed={processed} certified={certified} "
          f"uncertified<hmin={len(fails)} budget_exhausted={budget_exhausted} "
          f"time={time.time()-t0:.0f}s")
    if fails or budget_exhausted:
        print("=" * 72)
        print("VERDICT: FAIL -- det DF>0 NOT certified "
              + ("(uncertified boxes below hmin: try smaller --hmin / higher --dps)"
                 if fails else "(budget exhausted: raise --budget)"))
        sys.exit(1)
    print("=" * 72)
    print("VERDICT: PASS -- det DF(x) > 0 for all x in R^2, "
          "rigorously (validated intervals, exact witness.json input).")
    sys.exit(0)


if __name__ == "__main__":
    main()
