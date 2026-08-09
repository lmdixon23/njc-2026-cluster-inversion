# Asymptotic-polygon global inversion for planar saturating ridge networks and an injective neural-Jacobian separation witness

[![verify](https://github.com/lmdixon23/njc-2026-cluster-inversion/actions/workflows/verify.yml/badge.svg)](https://github.com/lmdixon23/njc-2026-cluster-inversion/actions/workflows/verify.yml)

This repository contains the manuscript source, exact data, and reproducibility
materials for Asymptotic polygon global inversion for planar saturating ridge
networks and an injective neural Jacobian separation witness.

Repository: https://github.com/lmdixon23/njc-2026-cluster-inversion

## Main result

The manuscript develops a global-inversion theorem for planar saturating ridge networks and resolves the canonical separation witness as an injective global diffeomorphism onto an explicit nonconvex octagon. Exact and validated parameter certificates show that the mechanism persists on an explicit open neighborhood.

## Repository contents

- `paper/` contains the LaTeX manuscript and its publication figure.
- `data/` and `config/` contain the exact witness and certificate inputs.
- `code/` contains maintained certificate, reconstruction, and figure scripts.
- `results/` contains the machine readable reports used by the manuscript.
- `verification/` maps claims to artifacts and documents the independent checks.
- `BLIND_CHECK.md` documents the separate reconstruction checks.

Local notes, exploratory work, build products, and superseded outputs belong in
the ignored `_local/` directory and are not part of the public repository.

## Maintained verification

Install the fully pinned Python environment and run the complete maintained gate:

```bash
python -m pip install -r requirements-lock.txt
bash run_all.sh
```

The expected final line is:

```text
VERDICT: ALL MAINTAINED CHECKS PASSED
```

The gate verifies complete repository integrity, exact witness and polygon
geometry, semantic consistency of every proof report and aggregate, the
structural analytic and winding package, independent Arb evidence, and focused
negative regressions for the previously identified false PASS paths.

## Related research

The separate [planar four-ridge counterexample repository](https://github.com/lmdixon23/njc-2026-planar-n4-winding) gives the exact positive-Jacobian noninjective counterexample at the next hidden width and the resulting sharp hidden-width threshold.

## Author

Logan M. Dixon · [research site](https://lmdixon23.github.io/) · [ORCID](https://orcid.org/0009-0001-0592-462X)

## Build

From the repository root:

```bash
python -m pip install -r requirements.txt
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/main.tex
```

The generated paper is `paper/main.pdf`.

## Component verification

```bash
python code/test_decimal_interval.py
python code/test_fail_closed.py
python code/verify_theorem_package.py
bash code/run_structural_checks.sh
```

The theorem package verifier recomputes report semantics, validates exact
component sets and witness, configuration, report, and leaf bindings, and
checks the shipped complete Arb replays. It does not rerun the full 93,000 leaf
and 96,517 leaf workflows during the fast gate.

## Full parameter certificate rerun

```bash
bash code/run_full_certificates.sh
```

This longer workflow regenerates the compact leaf partitions, tail reports,
complete Arb and Decimal replays, and semantically validated content addressed
combined reports.

The publication figure is generated directly from the canonical witness:

```bash
python code/generate_polygon_regimes.py
```

## Canonical witness

`data/witness.json` has SHA-256
`b4268c92e31066cfa3ea8466dc181ad3bf5876fc29b512493c5d37807c157a40`.

See `verification/CLAIM-TO-ARTIFACT-MAP.md` for the evidence map and
`RELEASE-NOTES.md` for the public release history.
