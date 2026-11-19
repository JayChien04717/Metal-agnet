---
name: sql-hfss-design
description: 'Router for the SQL HFSS repository (Design1-6 superconducting transmon chips, Qiskit Metal). Use for any layout, HFSS S21 / eigenmode, Q3D capacitance, Purcell filter and readout resonator, qubit Hamiltonian or CZ task in this repository. It points to the repository rules in AGENTS.md and skill/*.md. Do not use for the fluxonium FTF project or other repositories.'
---

# SQL HFSS design router

This skill only routes. The rules live in the repository, so read them there rather than relying on this file.

1. Read `AGENTS.md` at the repository root. It covers environments, run folders, provenance and confirmation before commit.
2. Read the root `README.md`, then the target `DesignN/README.md` and `DesignN/design.json`.
3. Then open the matching file in `skill/`:

| Task | Read |
|---|---|
| Parameters, results, provenance | `skill/parameter_provenance.md` |
| QComponents, layouts, DRC | `skill/skill.md`, `skill/qcomponents.md`, `skill/chip_drc.md` (plus `skill/launch_pad.md` for launch pads) |
| GDS export, figures | `skill/gds_export.md`, `skill/layout_figures.md` |
| Purcell filter, readout resonator, HFSS S21 fitting | `skill/purcell_resonator_design.md`, `note/purcell_filter_method.md` (QComponent `PurcellReadout`, circuit model `code/resonator/purcell_readout_model.py`; old pair parts are archived) |
| Two-qubit gate: CZ, ZZ off point, tunable coupler | `skill/two_qubit_gate_design.md`, `note/tunable_coupler_cz_method.md` |
| Capacitance matrix (Q3D), mutual inductance | `skill/cap_matrix_analysis.md`, `skill/mutual_inductance_analysis.md` |

4. Use the existing entry points under `DesignN/simulation/{general,qubit,resonator,flux}/` and `code/general/run_record.py`. Do not reconstruct scripts from another project.

Skills written for other projects do not apply here. If a skill names a different project path or interpreter (for example `ftf-simulate` and `FTF_project`), ignore it.
