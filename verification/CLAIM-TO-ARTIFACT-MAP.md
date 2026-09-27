# Claim to artifact map

| Paper item | Primary evidence | Separate check or failure condition |
|---|---|---|
| Negative and singular segment averages | `code/certify_avgdet_neg.py`, `data/witness.json` | `code/check_witness_core.py` and the continuity argument in the paper |
| Global determinant positivity for the witness | Analytic proof in `paper/main.tex`, `code/certify_analytic_determinant.py` | `code/check_analytic_det.py` and the validated compact and tail reports |
| Closed-width coordinate-dependence example | Analytic proof in `paper/main.tex` | `code/verify_closed_width_example.py` checks mixed minors, exact six-vertex hexagon geometry, area, simplicity, and the rational determinant bound |
| Jordan cluster set inversion theorem | Self contained proof in `paper/main.tex` | Every topological hypothesis is stated and checked in Appendix A |
| Asymptotic polygon construction | Proof in `paper/main.tex`, `code/asymptotic_polygon.py` | `code/polygon_indep.py` |
| Polygon classification and chamber structure | Exact derivation in `paper/main.tex`, `code/magnitude_chambers.py` | `code/verify_magnitude_chamber_refinement.py` |
| Witness octagon simplicity | `code/verify_witness_polygon.py` | `code/polygon_indep.py` and the exact nonadjacent pair test |
| Global injectivity and diffeomorphism conclusion | Theorem chain in `paper/main.tex` | Exact winding spectrum in `code/polygon_indep.py` |
| Degree equals winding | Proof in `paper/main.tex`, `code/polygon_winding.py` | Three exact arithmetic winding implementations |
| Open injective family | Proof in `paper/main.tex`, `code/certify_open_family.py` | Semantically validated and content addressed reports checked by `code/verify_theorem_package.py` |
| Common and anisotropic parameter boxes | Compact and tail reports under `results/` | Exact partition replay, complete separately implemented Arb leaf replay, and separate Decimal diagnostics |
| Fixed-radius tail grid | `results/point_tail_recomputed.json` | Complete separately implemented Arb replay in `results/point_tail_arb.json` |
| Release gate failure behavior | `code/report_validation.py`, `code/verify_sha256_manifest.py` | Adversarial regressions in `code/test_fail_closed.py` |

`verification/CLAIMS.md` records the public claim status. `bash run_all.sh` is
the maintained release gate. Full parameter certificate regeneration is kept
separate because it is substantially more expensive.
