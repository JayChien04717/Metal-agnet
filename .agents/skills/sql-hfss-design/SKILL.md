---
name: sql-hfss-design
description: 'Router for the SQL HFSS repository (Design1–Design6 superconducting qubit chips — fluxonium FTF, transmon, Purcell readout — Qiskit Metal). Use for any layout, HFSS S21 / eigenmode, Q3D capacitance, Purcell filter and readout resonator, qubit Hamiltonian or CZ task in this repository. It points to the repository rules in AGENTS.md and skill/*.md. Do not use for the separate FTF_project repository or other repositories.'
---

# SQL HFSS design router

This skill only routes. The rules live in the repository, so read them there rather than relying on this file.

1. Read `AGENTS.md` at the repository root. It covers the 3-stage gate rule (Stage 1 draft layout & figure review → Stage 2 EM simulation → Stage 3 fab GDS/e-beam export), environments, run folders, provenance and confirmation before commit. Never jump ahead across stages.
2. Read the root `README.md`, then the target `DesignN/README.md` and `DesignN/design.json`.
3. Then open the matching file in `skill/`:

| Task | Read |
|---|---|
| Parameters, run storage, provenance, confirmation | `skill/parameter_provenance.md` |
| QComponents, layouts, DRC, airbridges | `skill/skill.md`, `skill/layout_structure.md`, `skill/qcomponents.md`, `skill/chip_drc.md`, `skill/airbridge.md` |
| EM simulation environment, solver defaults, wire bonds | `skill/em_simulation.md` |
| Single-die GDS export (launch pads, cheesing, coupons), figures | `skill/gds_export.md`, `skill/layout_figures.md` |
| E-beam junction masks, EJ → resistance → junction area | `skill/ebeam_junction.md` |
| 4-inch wafer assembly & fab deliveries (laser writer, e-beam sites, markers) | `skill/fabrication_output.md` |
| Purcell filter, readout resonator, HFSS S21 fitting | `skill/purcell_resonator_design.md`, `note/purcell_filter_method.md` |
| Two-qubit gate: CZ, ZZ off point, tunable coupler | `skill/two_qubit_gate_design.md`, `note/tunable_coupler_cz_method.md` |
| Capacitance matrix (Q3D), mutual inductance | `skill/cap_matrix_analysis.md`, `skill/mutual_inductance_analysis.md` |

4. Use the existing entry points under `DesignN/simulation/{general,qubit,resonator,flux}/` and `code/general/run_record.py`. Do not reconstruct scripts from another project.

Skills written for other projects do not apply here. If a skill names a different project path or interpreter (for example `ftf-simulate` and `FTF_project`), ignore it.
