# Parameter provenance, run storage, and confirmation rules

Applies to new Designs, new runs, and organization of older results. Do not overwrite historical raw data; add provenance through sidecars without inventing missing evidence.

## 1. Roles and states

- **role:** `target`, `input`, `derived` (calculated/fitted), or `result`. Downstream inputs must retain upstream provenance.
- **status:** `draft`, `recorded` (stored does not mean validated), `validated` (supported by evidence), `unknown`, or `not_applicable`.
- **applicability:** `current`, `historical`, `needs_revalidation`, or `unknown`. Track validation status separately from geometry applicability.

## 2. Source boundaries

| Parameter | Source and requirements |
|---|---|
| EJ, EL | User settings or quantum-model solutions under explicit targets/constraints; record whether fixed or adjustable |
| Lumped Lj / Larray | EM inputs; preserve the conversion direction to/from EJ/EL |
| Capacitance matrix | Raw Q3D result; record net order, reference ground, units, and elimination |
| EC, charge coupling | Matrix derivation or EM circuit fitting; do not label as direct EM outputs |
| Linear modes, S/Z | HFSS; record boundaries, ports, mesh, frequency range, and convergence |
| Qubit f01, α, ZZ | Quantum model; record EJ/EL/EC, bias, and truncation |
| Fidelity, leakage | Dynamics calculation; record pulse, time step, algorithm, dissipation, and definitions |

- HFSS linear frequencies are not nonlinear qubit f01. If `scqubits` solves for EJ using HFSS-fitted EC, record it as `"scqubits target solve, dependent on HFSS-fitted EC"`, not as EJ obtained from HFSS.
- Never automatically retarget user-fixed EJ/EL. Mark unclear historical authorization `unknown`. If a transmon model has no independent EL, record `not_applicable` rather than `0`.

## 3. Per-parameter records

Required fields: `value`, `unit`, `role`, `source_kind`, `source_ref` (exact file and field/function), `method`, `depends_on`, `constraints`, `status`, `applicability`, `geometry_revision`, and `run_id`. Unknown values may be `null` with an explanation.

Explicitly label E/h in GHz. Distinguish f from ω and κ from κ/2π. Specify the Hamiltonian coefficient definition for coupling (charge vs. exchange). Φ/Φ0 is dimensionless.

## 4. Directory structure and per-run records

- **Shared tools:** `code/general/`, `code/qubit/`, `code/resonator/`, `code/flux/`.
- **Per-Design categories:** `layout/`, `simulation/{general,qubit,resonator,flux}/`, and matching `results/{general,qubit,resonator,flux}/`. Do not add `jobs/`, `support/`, `release/`, or `sweep/` layers.
- **New runs (`code/general/run_record.py`):** under a unique run ID, store `run.json`, `code/`, `raw/`, and `summaries/` / `figures/` / `logs/` (plus `inputs/` when input snapshots are needed). Record time/timezone, command, working directory, environment/versions, geometry revision, commit/dirty state, input hashes, upstream runs, execution status, and validation evidence.
- A commit hash alone does not identify a dirty workspace: save actual input snapshots.
- **Re-analysis:** create a separate run referencing `upstream` instead of overwriting adopted data.
- **Legacy data:** keep old mixed results in `general/legacy/` without rewriting original records (older migrated runs may still contain `manifest.json` / `analysis/`). Connect legacy entry points that write fixed filenames to a new run output directory before execution; a migration index is not historical provenance evidence.
- **Solver side files:** never write `hfss_eig_f_convergence.csv` or similar files to the repository root (`project_paths.configure()` and `run_record.create_run()` redirect them to `raw/`).

## 5. Adopting results

- Select `adopted_results` explicitly in `design.json` rather than picking the newest filename. README presents rounded summaries linked to exact source values.
- Mark earlier results `needs_revalidation` after changes to geometry, junction models, or boundaries. Restore `current` only with a rerun or concrete applicability check.

## 6. Final version confirmation before commit/push

- Before committing or pushing, prepare a reviewable summary of the final candidate, parameters, adopted runs, GDS/figures, and unvalidated items; obtain the user's explicit confirmation of the final version.
- Reuse explicit confirmation while version and parameters remain unchanged; any adjustment invalidates it and requires reconfirmation.
- Version confirmation does not itself authorize `git commit` / `git push`; perform Git operations only when requested by the user.
- README parameters and results must match only the final confirmed version. Keep older values in `notes/` or historical provenance indexes; mark missing new-version results as pending validation instead of substituting old data.
- Record `final_selection` (version, confirmation source, scope) and `release_confirmation` (whether the complete parameter snapshot is confirmed) in `design.json`. Selecting a candidate does not confirm every parameter, result, or fabrication revision.
