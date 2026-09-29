#!/usr/bin/env python3
"""Reject optimized execution of assertion-bearing verification modules."""
from pathlib import Path
import ast
import os
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class VerificationRuntimeTests(unittest.TestCase):
    def test_assertion_modules_reject_optimization(self):
        modules = [p for p in sorted((ROOT / "code").glob("*.py"))
                   if any(isinstance(n, ast.Assert)
                          for n in ast.walk(ast.parse(p.read_text(encoding="utf-8"))))]
        self.assertTrue(modules)
        for path in modules:
            for mode in ("-O", "-OO", "environment"):
                with self.subTest(module=path.name, mode=mode):
                    env = os.environ.copy()
                    env.pop("PYTHONOPTIMIZE", None)
                    env["PYTHONDONTWRITEBYTECODE"] = "1"
                    flags = [] if mode == "environment" else [mode]
                    if mode == "environment":
                        env["PYTHONOPTIMIZE"] = "1"
                    result = subprocess.run([sys.executable, "-B", *flags, str(path)],
                                            cwd=ROOT, env=env, capture_output=True,
                                            text=True, timeout=15)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("verification requires assertions", result.stderr)


if __name__ == "__main__":
    unittest.main()
