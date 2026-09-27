# Release notes

## Candidate v1.2.0 - 2026-09-28

- Corrected the closed-width coordinate-dependence example: the asymptotic hexagon now lists exactly its six attainable vertices, and the negative averaged-determinant claim is proved by an explicit rational bound rather than an informal asymptotic estimate.
- Rewrote the generic self-crossing argument in terms of signed multiplicity $\varepsilon\,\mathrm{wind}$, eliminating an orientation-sign ambiguity while preserving the criterion.
- Added a fail-closed exact regression for the closed-width example and connected it to the maintained theorem-package gate.
- Tightened verification provenance language to distinguish separate internal implementations from external review, refreshed bibliographic metadata, and removed a nonessential higher-dimensional aside.
- Regenerated the fixed-radius Arb tail replay with wall-clock timing removed from the persisted JSON, making that canonical certificate deterministic under exact rerun.
- Refreshed manuscript date and prospective `v1.2.0` release references for the publication candidate; the tag must be minted only after the exact reviewed bytes are finalized.

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
