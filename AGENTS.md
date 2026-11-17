# Project AI Working Rules

## Task entry points

- Read the root README, the target Design README, and design.json when available.
- For parameters or results, read skill/parameter_provenance.md.
- For QComponents or layouts, read skill/skill.md, skill/qcomponents.md, and skill/chip_drc.md; for launch pads, also read skill/launch_pad.md.
- For GDS or figures, also read skill/gds_export.md and skill/layout_figures.md.
- For Purcell filter / readout resonator channels, read skill/purcell_resonator_design.md and note/purcell_filter_method.md.
- For two-qubit gates (CZ, ZZ, tunable coupler, coupler flux), read skill/two_qubit_gate_design.md and note/tunable_coupler_cz_method.md.
- For capacitance or mutual-inductance analysis, read skill/cap_matrix_analysis.md or skill/mutual_inductance_analysis.md, respectively.
- Copy and complete templates/Design/ for new designs; the template is not a finished design.

## Skill routing

- The operating rules for this repository are the files in `skill/` (entry points above). The registered router skill is `.agents/skills/sql-hfss-design/SKILL.md`, which points back to them.
- Global or user-level skills written for other projects do not apply here, even when the task mentions HFSS, S21, Q3D, resonators or capacitance. One example is `ftf-simulate`, which is for the fluxonium FTF project and now lives in that project's folder.
- If a selected skill names another project path, environment or design, stop and use this repository's `skill/` files instead.

## Organization

- Keep shared QComponents in pattern design/, new chip assemblies in layout/, and design-specific components in layout/components/.
- Keep operating rules in skill/, general methods in note/, and design-specific values, exceptions, and history within each Design.
- All Designs are organized by purpose. code/general/project_paths.py and the migration index resolve cross-design imports and historical result paths. New code must use the new locations directly.
- Preserve unrelated uncommitted changes; do not overwrite or revert them.

## Parameters and evidence

- Do not automatically retarget user-fixed EJ/EL. Record any gap between the constraints and the requested targets.
- EJ/EL are user or quantum-model design inputs; HFSS lumped L is also an input. Do not describe them as quantities independently extracted by HFSS.
- Distinguish target, input, derived, and result. Record provenance, units, assumptions, and validation status for each parameter.
- Select adopted runs/files explicitly in design.json rather than choosing the most recent timestamp.
- Mark missing sources, authorization, commits, hashes, or convergence evidence as unknown; do not invent them.
- After geometry changes, mark earlier results historical or needs_revalidation.

## Execution and delivery

- Use the Design-specified qmetal environment for layout/EM work and scqenv or another documented environment for quantum analysis.
- Preserve existing solver defaults. Pass counts are execution settings, not convergence evidence; record actual checks or missing evidence.
- Preserve immutable input snapshots, raw data, analysis, environment details, and status for new runs. Save re-analysis separately and reference its parent run.
- Label simulation and fabrication geometries separately. GDS and figures must correspond to the same revision.
- Keep current summaries, sources, reproduction entry points, and limitations in README; keep detailed history in notes.
- Check links and JSON after documentation changes. Explicitly state when simulations have not been rerun.

This document and the provenance rules supersede the old conventions of using one environment for every script and reading the newest parameter file. Other existing geometry and fabrication rules remain applicable.


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

## Documentation language

Write all skills, notes, and README files in English, including new documents and future updates. Preserve technical identifiers, source paths, numerical values, units, and evidence status when translating. Historical raw data and immutable run snapshots retain their original content.
