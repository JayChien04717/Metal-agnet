# Project AI Working Rules

## Task entry points

- Read the root README, the target Design README, and `design.json` when available (including any Design-specific rules or solver overrides).
- For parameters, run storage, or final confirmation, read `skill/parameter_provenance.md`.
- For QComponents or layouts, read `skill/skill.md`, `skill/layout_structure.md`, `skill/qcomponents.md`, and `skill/chip_drc.md`.
- For EM simulation environments, solver defaults (Q3D 8 passes; HFSS eigenmode 6 passes; HFSS driven 4 passes with explicit mesh length), and wire bonds, read `skill/em_simulation.md`.
- For single-die GDS (including launch pads, cheesing, test coupons) or layout figures, read `skill/gds_export.md` and `skill/layout_figures.md`.
- For 4-inch wafer assembly and deliveries sent to fabrication (laser writer, NTHU / TSRI / NCKU e-beam), read `skill/fabrication_output.md`.
- For Purcell filter / readout resonator channels, read `skill/purcell_resonator_design.md` and `note/purcell_filter_method.md`; build new channels with `PurcellReadout` (`pattern design/resonator/purcell_readout.py`; circuit model `code/resonator/purcell_readout_model.py`).
- For e-beam junction masks (small JJ, JJ array, SQUID) or EJ → resistance → area → mask width, read `skill/ebeam_junction.md`; use `FluxoniumJunctionDolan` / `TransmonJunctionDolan` in `pattern design/ebeam_junction/` and `code/qubit/jj_resistance_model.py`.
- For two-qubit gates (CZ, ZZ, tunable coupler, coupler flux), read `skill/two_qubit_gate_design.md` and `note/tunable_coupler_cz_method.md`.
- For capacitance or mutual-inductance analysis, read `skill/cap_matrix_analysis.md` or `skill/mutual_inductance_analysis.md`, respectively.
- Copy and complete `templates/Design/` for new designs; the template is not a finished design.

## Skill routing

- The operating rules for this repository are the files in `skill/` (entry points above). The registered router skill is `.agents/skills/sql-hfss-design/SKILL.md`, which points back to them.
- Global or user-level skills written for other projects do not apply here, even when the task mentions HFSS, S21, Q3D, resonators or capacitance. One example is `ftf-simulate`, which is for the fluxonium FTF project and now lives in that project's folder.
- If a selected skill names another project path, environment or design, stop and use this repository's `skill/` files instead.

## Organization

- Keep shared QComponents in `pattern design/`, new chip assemblies in `layout/`, and design-specific components in `layout/components/`.
- Keep `skill/` (general operating rules) and `note/` (general physics/methods) strictly **design-agnostic**: never hardcode specific Design names, numbers, or Design-specific exceptions there. Put Design-specific values, solver overrides (e.g., reduced passes for large multi-channel chips), script entry points, and history inside each `Design<N>/` folder (`README.md`, `notes/`, `design.json`).
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

- **Follow the 3-stage gate workflow (`skill/skill.md` §2) — never jump ahead:**
  1. **Stage 1 (Draft layout & figure review):** Build only simulation `QComponent`s and `layout/build.py`, verify basic layout DRC (pins outward, filleted lengths on target, no crossings, component-level 4 µm metal width), and render `fig/chip.png` + part figures from `*_sim.gds` ONLY. **Stop and wait for user review.** Do not create `*_fab.py` variants, place `FluxoniumJunctionDolan` / `TransmonJunctionDolan` e-beam layers, export `*_fab.gds` / test coupons / cheesing, or run EM simulations unless explicitly asked.
  2. **Stage 2 (EM simulation & tuning):** Run Q3D / HFSS only after the user confirms the layout (or explicitly asks to simulate), iterating geometry until parameters are confirmed.
  3. **Stage 3 (Fab GDS & wafer delivery):** Only after simulation parameters are confirmed (or the user explicitly requests fab output), ask the `skill/gds_export.md` §2.1 checklist questions, build `*_fab` / e-beam layers, export `*_fab.gds`, and run full mask checks.
- Use the Design-specified qmetal environment for layout/EM work and scqenv or another documented environment for quantum analysis.
- Preserve existing solver defaults. Pass counts are execution settings, not convergence evidence; record actual checks or missing evidence.
- Preserve immutable input snapshots, raw data, analysis, environment details, and status for new runs. Save re-analysis separately and reference its parent run.
- Label simulation and fabrication geometries separately. GDS and figures must correspond to the same revision.
- Keep current summaries, sources, reproduction entry points, and limitations in README; keep detailed history in notes.
- Check links and JSON after documentation changes. Explicitly state when simulations have not been rerun.

This document and the provenance rules supersede the old conventions of using one environment for every script and reading the newest parameter file. Other existing geometry and fabrication rules remain applicable.


## Final version confirmation and storage conventions

Full rules live in `skill/parameter_provenance.md` (§4 and §6):
- **Before commit/push:** prepare a reviewable summary (final candidate, parameters, adopted runs, GDS/figures, unvalidated items) and obtain explicit user confirmation; record `final_selection` and `release_confirmation` in `design.json`. Confirmation does not itself authorize Git operations.
- **Directory and run layout:** shared tools in `code/{general,qubit,resonator,flux}/`; per-Design `layout/`, `simulation/{general,qubit,resonator,flux}/`, and matching `results/`. New runs store `run.json`, `code/`, `raw/`, and `summaries/` / `figures/` / `logs/` (plus `inputs/` when needed); never overwrite historical data or write solver side files to the repository root.
- **What GitHub tracks** (`skill/parameter_provenance.md` §7): code (`layout/`, `simulation/`, shared tools), `design.json` / `parameters/`, and only the records (`run.json`, `summaries/`, `figures/`) of runs `design.json` references. GDS files and all other run data stay local; after changing adopted runs, run `code/general/sync_tracked_results.py --apply` before committing.

## Documentation language

Write all skills, notes, and README files in English, including new documents and future updates. Preserve technical identifiers, source paths, numerical values, units, and evidence status when translating. Historical raw data and immutable run snapshots retain their original content.
