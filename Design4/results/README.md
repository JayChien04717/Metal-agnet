# Design4 simulation results

[File index](index.csv): category, original path, new location, and SHA-256. Original manifests and source snapshots are retained; a migration index is not evidence of the original execution version.

- [resonator/](resonator/)

qubit: capacitance/spectra/gates; resonator: readout/Purcell; flux: mutual inductance/pulses; general: cross-purpose summaries and checks.

legacy preserves historical data. Older names such as latest/final do not imply current user adoption. New executions use independent run directories and never overwrite legacy. See the Design's design.json for adopted results.

2026-11-16 cleanup: runs not cited by this Design's README.md, design.json or parameters/, or by 20261115_latest/README.md, were condensed into a single `parameter_report.json` (uncited loose legacy files into `legacy_parameter_report.json`) holding the run metadata, parameters, small result tables and a file list. Full raw data, code snapshots and figures of those runs remain in git commit `136d032fa62a`; local AEDT projects/solutions were deleted and are not recoverable. Cited runs are kept in full except their AEDT files.
