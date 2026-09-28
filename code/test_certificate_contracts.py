#!/usr/bin/env python3
"""Bounded regressions for certificate scope, evidence and resume contracts."""
import json
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from fractions import Fraction as F
from pathlib import Path

from mpmath import iv, mp
import certify_parameter_radius as cp
import report_validation as rv
from certify_tail_parameter_transfer import worst_positive_pair

ROOT = Path(__file__).resolve().parents[1]


class CertificateContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scratch = ROOT / '_local'
        scratch.mkdir(exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(prefix='contracts-', dir=scratch)
        cls.temp = Path(cls.temporary.name)
        cls.witness = ROOT / 'data/witness.json'
        cls.A, cls.B, cls.c, cls.sha = cp.load_witness(cls.witness)
        cls.cfg = rv.load_json(ROOT / 'config/witness_common.json')
        mp.dps = iv.dps = 80

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_existing_canonical_reports(self):
        for name, config in [('witness_common_5e-6', 'witness_common'), ('witness_anisotropic', 'witness_anisotropic')]:
            rv.validate_combined(ROOT / 'results' / name, ROOT / 'config' / (config + '.json'))

    def test_width_five_is_rejected(self):
        A = [self.A[0] + [F(0)], self.A[1] + [F(-1000)]]
        B = self.B + [[F(0), F(1)]]
        with self.assertRaises(ValueError):
            AI, BI, CI = cp.parameter_intervals(A, B, self.c + [F(0)], F(0))
            cp.interior_certify(1, BI, CI, cp.minor_products_iv(AI, BI), budget=100)

    def test_invalid_grid_rejected(self):
        AI, BI, CI = cp.parameter_intervals(self.A, self.B, self.c, F(0))
        M = cp.minor_products_iv(AI, BI)
        pos, neg, _ = cp.sign_structure(M)
        with self.assertRaises(ValueError):
            cp.tail_certify(582, BI, CI, M, pos, neg, Mgrid=-1)

    def test_skip_everything_rejected(self):
        with self.assertRaises(ValueError):
            cp.certify('0', self.witness, 582, 200000, F('.008'), 100, 0,
                       skip_tail=True, skip_interior=True)

    def test_same_context_bad_resume_rejected(self):
        c = [F(1000), F(1000), F(0), F(0)]
        AI, BI, CI = cp.parameter_intervals(self.A, self.B, c, F(0))
        M = cp.minor_products_iv(AI, BI)
        path = '10101'
        bad_path = self.temp / 'resume-witness.json'
        bad_path.write_text(json.dumps({'A': [[str(x) for x in row] for row in self.A],
                                       'B': [[str(x) for x in row] for row in self.B],
                                       'c': list(map(str, c))}), encoding='utf-8')
        context = cp.certificate_context('interior', rv.file_sha256(bad_path), F(0), F(0), F(0),
                                         582, 200000, F('.008'), 40000000, 80)
        state = {'schema': 'parameter-radius-interior-state-v2', 'context': context,
                 'context_sha256': cp.canonical_hash(context), 'processed': 10, 'certified': 5,
                 'paths': ['0', '11', '100', '1011', '10100'], 'worst': '1',
                 'stack': [{'path': path, 'box': list(map(str, cp.path_box_exact(path, F(582))))}]}
        state_path = self.temp / 'resume.json'
        state_path.write_text(json.dumps(state), encoding='utf-8')
        with self.assertRaises(ValueError):
            cp.interior_certify(582, BI, CI, M, state_file=str(state_path), state_context=context)

    def test_negative_K_does_not_replace_corrected_slope(self):
        # Raw g=1 and K=-10 make R*g-K positive, but L*pi/grid>1
        # leaves a negative corrected slope, invalid for r -> infinity.
        def rho(i, cosine, sine):
            return iv.mpf(1) if i == 2 else iv.mpf(0)
        setup = ((2, 3), rho, {(0, 1): iv.mpf(-10)}, {(0, 1): mp.mpf(1)})
        with patch.object(cp, 'tail_setup', return_value=setup):
            ok, _, _ = cp.tail_certify(1, [[0, 0]] * 4, [0] * 4, {}, [(0, 1)], [(2, 3)], Mgrid=1)
        self.assertFalse(ok)

    def test_tail_worst_pair_uses_exact_ratios(self):
        half = F(1, 2)
        larger = half + F(1, 10**100)
        self.assertEqual(float(half), float(larger))
        minors = {(0, 1): F(1), (0, 3): F(1), (2, 3): F(-1)}
        errors = {(0, 1): half, (0, 3): larger, (2, 3): F(0)}
        self.assertEqual(worst_positive_pair(minors, errors, (2, 3)), (0, 3))

    def test_strict_bound_parsing(self):
        self.assertEqual(rv.first_decimal('[1 +/- 2]', 'ball'), F(-1))
        self.assertEqual(rv.first_decimal('[.1, .2]', 'interval'), F(1, 10))
        for value in ('nan 1', 'garbage123', 'Infinity', '[2, 1]', '[1 +/- -2]'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                rv.first_decimal(value, 'invalid')

    def test_base_margin_threshold(self):
        point = rv.load_json(ROOT / 'results/point_tail_recomputed.json')
        point['tail']['worst_margin'] = '1e-100'
        with self.assertRaises(ValueError):
            rv.validate_point_tail(point, self.sha, 582, 200000)

    def test_required_geometry_predicates(self):
        with self.assertRaises(ValueError):
            rv.validate_geometry({'verdict': 'PASS', 'witness_sha256': self.sha,
                                  'checks': {'unrelated': True}, 'segment_failures': []}, self.sha)

    def test_geometry_radius_must_cover_box(self):
        geometry = rv.load_json(ROOT / 'results/witness_common_5e-6/geometry.json')
        with self.assertRaises(ValueError):
            rv.validate_geometry(geometry, self.sha, dict(self.cfg, delta_A='1'))

    def test_missing_volume_rejected(self):
        report = rv.load_json(ROOT / 'results/witness_common_5e-6/compact_independent.json')
        report.pop('volume')
        with self.assertRaises(ValueError):
            rv.validate_compact_independent(report, self.cfg, {'count': 93000})

    def test_negative_arb_lower_rejected(self):
        folder = ROOT / 'results/witness_common_5e-6'
        report = rv.load_json(folder / 'compact_arb.json')
        report['minimum_leaf_lower'] = '[1 +/- 2]'
        with self.assertRaises(ValueError):
            rv.validate_compact_arb(report, self.cfg, ROOT / 'config/witness_common.json', self.sha,
                                    folder / 'compact.json', folder / 'compact_leaf_paths.json', {'count': 93000})

    def test_arb_base_threshold_rejected(self):
        report = rv.load_json(ROOT / 'results/point_tail_arb.json')
        report['worst_margin_lower'] = '[1e-100 +/- 1e-120]'
        path = self.temp / 'small-arb-base.json'
        path.write_text(json.dumps(report), encoding='utf-8')
        with self.assertRaises(ValueError):
            rv.validate_point_tail_arb(path, ROOT / 'results/point_tail_recomputed.json', self.witness, 582, 200000)

    def test_legacy_diagnostic_cannot_be_rebound(self):
        report = rv.load_json(ROOT / 'results/witness_common_5e-6/tail_independent.json')
        report['unverified_extra'] = True
        path = self.temp / 'rebound.json'
        path.write_text(json.dumps(report), encoding='utf-8')
        with self.assertRaises(ValueError):
            rv.validate_component_context('tail_independent', report, path, self.cfg,
                ROOT / 'config/witness_common.json', self.sha,
                ROOT / 'results/witness_common_5e-6/compact_leaf_paths.json',
                rv.validated_tail_base(self.witness, 582, 200000))

    def test_new_transfer_and_decimal_contexts(self):
        for config in ('witness_common', 'witness_anisotropic'):
            cfg = rv.load_json(ROOT / 'config' / (config + '.json'))
            da, db, dc = rv.config_deltas(cfg)
            # CLI accepts decimal radii, so retain the exact source tokens.
            common = cfg.get('delta', '0')
            args = ['--witness', str(self.witness), '--delta-a', cfg.get('delta_A', common),
                    '--delta-b', cfg.get('delta_B', common), '--delta-c', cfg.get('delta_c', common)]
            for script, name, option in [('certify_tail_parameter_transfer.py', 'tail', '--out'),
                                         ('independent_tail_transfer.py', 'tail_independent', '--json-out')]:
                path = self.temp / (config + '-' + name + '.json')
                run = subprocess.run([sys.executable, '-B', str(ROOT / 'code' / script), *args, option, str(path)],
                                     capture_output=True, text=True, cwd=ROOT)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                report = rv.load_json(path)
                if name == 'tail':
                    rv.validate_tail(report, cfg, self.sha, rv.abc_payload_sha256(self.witness))
                else:
                    rv.validate_tail_independent(report, cfg)
                    self.assertEqual(report['evidence_role'], 'diagnostic')
                    self.assertNotIn('residual_lower', report)
                rv.validate_component_context(name, report, path, cfg, ROOT / 'config' / (config + '.json'),
                    self.sha, ROOT / 'results' / cfg['name'] / 'compact_leaf_paths.json',
                    rv.validated_tail_base(self.witness, 582, 200000))

    def test_exact_payload_decimal_token(self):
        import hashlib
        value = '0.100000000000000000000000000001'
        path = self.temp / 'precise.json'
        path.write_text('{"A":[[' + value + ']],"B":[[0]],"c":[0]}', encoding='utf-8')
        payload = json.dumps({'A': [[value]], 'B': [['0']], 'c': ['0']}, sort_keys=True, separators=(',', ':')).encode()
        self.assertEqual(rv.abc_payload_sha256(path), hashlib.sha256(payload).hexdigest())

    def test_wrong_transfer_radius_rejected(self):
        run = subprocess.run([sys.executable, '-B', str(ROOT / 'code/certify_tail_parameter_transfer.py'),
                              '--witness', str(self.witness), '--delta', '0', '--R', '1', '--grid', '1'],
                             capture_output=True, text=True, cwd=ROOT)
        self.assertNotEqual(run.returncode, 0, run.stdout)

    def test_changed_bias_transfer_rejected(self):
        bad = json.loads(self.witness.read_text())
        bad['c'] = [1000, 1000, 0, 0]
        path = self.temp / 'bad-witness.json'
        path.write_text(json.dumps(bad), encoding='utf-8')
        run = subprocess.run([sys.executable, '-B', str(ROOT / 'code/certify_tail_parameter_transfer.py'),
                              '--witness', str(path), '--delta', '0'], capture_output=True, text=True, cwd=ROOT)
        self.assertNotEqual(run.returncode, 0, run.stdout)


if __name__ == '__main__':
    unittest.main(verbosity=2)
