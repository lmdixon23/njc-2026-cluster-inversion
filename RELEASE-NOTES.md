# Release notes

## v1.1.0 - 2026-07-27

- Added a section proving the closed widths in graph coordinates and a proposition showing the averaged-Jacobian criterion is coordinate dependent.
- Reported the optimizer clamp and bridge loss for the analytic certificate; corrected the anisotropic box description; stated the independence of the Arb replay; bounded the novelty claim.
- Fixed cross-platform reproducibility: LF checkout and LF artifact writers, timing-free persisted reports, ASCII gate output, and dependency-aware interpreter selection in run_all.sh.

## v1.0.0 — 2026-07-17

- Provides the self contained Paper II manuscript, including the separation
  result needed from Paper I and the global injectivity resolution.
- Includes the exact witness, analytic and interval certificates, asymptotic
  polygon checks, winding computations, and parameter box reports.
- Includes a fully pinned Python environment, a one command verification gate,
  a GitHub Actions workflow, a public claim map, and a complete checksum
  manifest.
- Semantically validates every compact, tail, geometry, replay, and aggregate
  report rather than trusting top-level verdict fields.
- Adds public complete Arb replays for all compact leaves and all fixed-radius
  angular directions, bound to the exact witness and source reports.
- Binds compact resume states to their full context and exact prefix-free root
  coverage, disables unsafe tail resume, and returns status 2 for paused work.
- Adds focused negative regressions for omitted analytic premises, forged
  checkpoints, contradictory reports, missing aggregate components, and
  incomplete manifests.
- Generates the publication figure from canonical exact data with embedded
  TrueType fonts and a deterministic output path.
- Uses repository relative paths and a portable PDFLaTeX build.
