#!/usr/bin/env python3
"""The symbolic producer must reject predicates outside its proof contract."""
import unittest
import sympy as sp
from magnitude_chambers import primitive_coefficients


class MagnitudePredicateTests(unittest.TestCase):
    def test_rejects_nonlinear_or_inhomogeneous_predicates(self):
        x, y = sp.symbols("x y")
        for expression in (x*x, x*y, x + 1, sp.Integer(1)):
            with self.subTest(expression=expression):
                with self.assertRaises(ValueError):
                    primitive_coefficients(expression, (x, y))

    def test_exact_linear_normalization(self):
        x, y = sp.symbols("x y")
        self.assertEqual(primitive_coefficients(-x/2 + y/3, (x, y)), (3, -2))
        self.assertIsNone(primitive_coefficients(sp.Integer(0), (x, y)))


if __name__ == "__main__":
    unittest.main()
