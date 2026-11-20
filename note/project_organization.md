# Project organization

General directory and storage layout for all Designs in this repository. Operational run-recording and confirmation rules live in [skill/parameter_provenance.md](../skill/parameter_provenance.md).

- **Functional categories:** code and results are grouped into `general`, `qubit`, `resonator`, and `flux`; chip layouts (`layout/`) and shared `QComponent`s (`pattern design/`) are kept separate.
- **Shared vs. Design-specific:**
  - Reusable utilities and models belong in root `code/`.
  - Shared `QComponent`s belong in root `pattern design/`.
  - Scripts with Design constants, Design-specific component forks (`layout/components/`), and Design-specific solver overrides stay inside each `Design<N>/` folder.
- **Legacy and immutable records:**
  - Pre-reorganization mixed batches remain under `results/<category>/legacy/` without rewriting original manifests or snapshots; `results/index.csv` maps historical paths and hashes.
  - Local dated delivery folders at the repository root and under `Fabrication output/<YYYYMMDD>/` are git-ignored.
  - Reorganization record and hash validation: [`code/general/migration_record.json`](../code/general/migration_record.json) and [`code/general/migration_validation.json`](../code/general/migration_validation.json).
