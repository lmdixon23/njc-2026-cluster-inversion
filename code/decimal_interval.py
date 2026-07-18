#!/usr/bin/env python3
"""Small fail-closed Decimal interval arithmetic module.

This implementation is intentionally independent of mpmath.iv.  Transcendental
endpoints are evaluated at high Decimal precision and expanded by an explicit
absolute/relative pad.  It is a corroborating implementation, not a replacement
for a formally verified interval library.
"""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

PREC = 70
PAD = Decimal('1e-1000')


def D(x) -> Decimal:
    return x if isinstance(x, Decimal) else Decimal(str(x))


def _down(fn):
    with localcontext() as ctx:
        ctx.prec = PREC
        ctx.rounding = ROUND_FLOOR
        v = fn()
        return v - PAD * (Decimal(1) + abs(v))


def _up(fn):
    with localcontext() as ctx:
        ctx.prec = PREC
        ctx.rounding = ROUND_CEILING
        v = fn()
        return v + PAD * (Decimal(1) + abs(v))


@dataclass(frozen=True)
class DI:
    lo: Decimal
    hi: Decimal

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('invalid interval')

    @staticmethod
    def point(x) -> 'DI':
        x = D(x)
        return DI(x, x)

    def __add__(self, other):
        o = as_di(other)
        return DI(_down(lambda: self.lo + o.lo), _up(lambda: self.hi + o.hi))

    __radd__ = __add__

    def __neg__(self):
        return DI(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_di(other))

    def __rsub__(self, other):
        return as_di(other) - self

    def __mul__(self, other):
        o = as_di(other)
        products = [(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        lows = [_down(lambda a=a, b=b: a * b) for a, b in products]
        highs = [_up(lambda a=a, b=b: a * b) for a, b in products]
        return DI(min(lows), max(highs))

    __rmul__ = __mul__

    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval contains zero')
        lo1 = _down(lambda: Decimal(1) / self.lo)
        lo2 = _down(lambda: Decimal(1) / self.hi)
        hi1 = _up(lambda: Decimal(1) / self.lo)
        hi2 = _up(lambda: Decimal(1) / self.hi)
        return DI(min(lo1, lo2), max(hi1, hi2))

    def __truediv__(self, other):
        return self * as_di(other).reciprocal()

    def exp(self):
        return DI(_down(lambda: self.lo.exp()), _up(lambda: self.hi.exp()))

    def abs(self):
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return DI(Decimal(0), max(-self.lo, self.hi) + PAD)


def as_di(x) -> DI:
    return x if isinstance(x, DI) else DI.point(x)


def sigma_prime(t: DI) -> DI:
    den = DI.point(2) + t.exp() + (-t).exp()
    return den.reciprocal()
