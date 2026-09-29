#!/usr/bin/env python3
"""Small regression tests for the custom Decimal interval kernel."""

# Assertions below are part of verification, not optional diagnostics.
import sys as _verification_sys
if _verification_sys.flags.optimize:
    raise SystemExit("FAIL: verification requires assertions; omit -O/-OO and unset PYTHONOPTIMIZE.")

from decimal import Decimal, localcontext
from decimal_interval import DI


def main() -> int:
    x = Decimal('1.0000000000000000000000000001')
    product = DI.point(x) * DI.point(x)
    with localcontext() as ctx:
        ctx.prec = 120
        exact_product = x * x
    assert product.lo <= exact_product <= product.hi

    reciprocal = DI.point(Decimal(3)).reciprocal()
    with localcontext() as ctx:
        ctx.prec = 120
        exact_reciprocal = Decimal(1) / Decimal(3)
    assert reciprocal.lo <= exact_reciprocal <= reciprocal.hi

    mixed = DI(Decimal('-2.1'), Decimal('3.7')) * DI(Decimal('-5.2'), Decimal('0.4'))
    exact_corners = [
        Decimal('-2.1') * Decimal('-5.2'),
        Decimal('-2.1') * Decimal('0.4'),
        Decimal('3.7') * Decimal('-5.2'),
        Decimal('3.7') * Decimal('0.4'),
    ]
    assert mixed.lo <= min(exact_corners)
    assert mixed.hi >= max(exact_corners)
    print('[PASS] Decimal interval regression tests')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
