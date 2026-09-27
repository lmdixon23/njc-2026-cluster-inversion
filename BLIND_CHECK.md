# Separate reconstruction checks

The following scripts use implementations methodologically separate from the main certificate
paths:

All of these are internal checks. Their methodological separation concerns code path, arithmetic, or proof reconstruction; it does not mean external human review. Legacy filenames containing `independent` use that word only in the implementation-separation sense.

```bash
python code/check_analytic_det.py
python code/check_witness_core.py
python code/polygon_indep.py
python code/replay_leaf_paths_arb.py config/witness_common.json
python code/replay_leaf_paths_arb.py config/witness_anisotropic.json
python code/replay_point_tail_arb.py
```

They reconstruct the analytic determinant margin, witness mixed minors and
averaged determinant, exact polygon winding spectrum, every exported compact
leaf, and every fixed-radius angular tail direction from the canonical input
data. Each script exits nonzero when a release-critical predicate fails.

These checks are corroborating implementations. The paper proofs and validated
interval certificates retain the evidence roles stated in
`verification/CLAIM-TO-ARTIFACT-MAP.md`.
