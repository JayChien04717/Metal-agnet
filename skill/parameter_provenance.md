# Parameter provenance and version rules

Applies to new Designs, new runs, and organization of older results. Do not overwrite historical raw data. Add provenance through sidecars without inventing missing evidence.

## Roles and states

- role: target, input, derived (calculated/fitted), or result. Downstream inputs must retain upstream provenance.
- status: draft, recorded (stored does not mean validated), validated (supported by evidence), unknown, or not_applicable.
- applicability: current, historical, needs_revalidation, or unknown. Track validation status separately from geometry applicability.

## Source boundaries

| Parameter | Source and requirements |
|---|---|
| EJ, EL | User settings or quantum-model solutions under explicit targets/constraints; record whether fixed or adjustable |
| Lumped Lj/Larray | EM inputs; preserve the conversion direction to/from EJ/EL |
| Capacitance matrix | Raw Q3D result; record net order, reference ground, units, and elimination |
| EC, charge coupling | Matrix derivation or EM circuit fitting; do not label as direct EM outputs |
| Linear modes, S/Z | HFSS; record boundaries, ports, mesh, frequency range, and convergence |
| Qubit f01, alpha, ZZ | Quantum model; record EJ/EL/EC, bias, and truncation |
| Fidelity, leakage | Dynamics calculation; record pulse, time step, algorithm, dissipation, and definitions |

HFSS linear frequencies are not nonlinear qubit f01. If scqubits solves for EJ using HFSS-fitted EC, record this as "scqubits target solve, dependent on HFSS-fitted EC," not as EJ obtained from HFSS.

Never automatically retarget user-fixed EJ/EL. Mark unclear historical authorization unknown. If a transmon model has no independent EL, record not_applicable rather than presenting 0 as a result.

## Per-parameter records

Required fields: value, unit, role, source_kind, source_ref (exact file and field/function), method, depends_on, constraints, status, applicability, geometry_revision, and run_id. Unknown values may be null with an explanation.

Explicitly label E/h in GHz. Distinguish f from omega and kappa from kappa/2pi. Specify the Hamiltonian coefficient definition for coupling; do not confuse charge coefficients with exchange coupling. Phi/Phi0 is dimensionless.

## Per-run records

Under a unique run ID, retain manifest.json, inputs/, raw/, and analysis/. Record time/timezone, command, working directory, environment/versions, geometry revision, commit/dirty state, input hashes, upstream runs, execution status, and validation evidence.

A commit alone does not identify a dirty workspace; save actual input snapshots. Reanalysis requires a new analysis ID referencing the original run and must not overwrite adopted data. Leave unknown template fields null; complete them before a new run and explicitly document irrecoverable gaps in older runs.

## Adopting results

Select adopted_results in design.json rather than choosing the newest filename. README provides rounded summaries linked to exact source values.

Mark results needs_revalidation after changes to geometry, junction models, or boundaries. Restore current only with a rerun or concrete applicability checks. An older report's claim of no impact is insufficient evidence on its own.


## Final version confirmation before commit/push

- Before committing or pushing, prepare a reviewable summary of the final candidate, parameters, adopted runs, GDS/figures, and unvalidated items; then obtain the user's confirmation of the final version.
- Reuse explicit confirmation while the version and parameters remain unchanged. Subsequent version or parameter adjustments invalidate that confirmation and require reconfirmation before submission.
- Version confirmation does not itself authorize a commit/push; perform Git operations as requested by the user.
- README parameters and results must correspond only to the final confirmed version. Keep older values in notes or historical provenance indexes; mark missing new-version results as pending validation instead of substituting old data.
- Record final_selection (version, confirmation source, scope) and release_confirmation (whether the complete parameter snapshot is confirmed) in design.json. Selecting a candidate does not confirm every parameter, result, or fabrication revision.


## Current directory and storage conventions (2026-11-15)

Follow the user's simplified categories: shared tools go in code/general or code/qubit, with resonator/flux added when needed. Each Design uses layout/ and simulation/{general,qubit,resonator,flux}/, with matching result categories. Do not add jobs, support, release, or sweep layers.

New runs use run.json, code/, raw/, and summaries/figures/logs. Record settings and provenance in run.json; add inputs/ when input snapshots are needed. Create a separate run for re-analysis and record upstream. Keep old mixed results in general/legacy without rewriting original records; a migration index is not historical provenance evidence.

Some legacy entry points still write fixed filenames. Connect them to a new run output directory before execution; do not overwrite historical results. Path compatibility does not mean every legacy entry point is automatically integrated with the run recorder.
