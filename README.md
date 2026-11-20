# Superconducting Quantum Chip Layout & Simulation Project

Planar superconducting qubit chip design pipeline: build layouts with **Qiskit Metal**, extract electromagnetic parameters with **Ansys Q3D / HFSS**, model Hamiltonians and two-qubit gates with **scqubits / QuTiP**, and export verified GDS masks and 4-inch wafer deliveries.

## Getting started

1. Read the [AI working rules](AGENTS.md) and [parameter provenance rules](skill/parameter_provenance.md).
2. Open the target `Design*/` directory and check its `README.md` and `design.json` for its current revision, adopted runs, execution environment, and any Design-specific overrides.
3. To start a new chip design, copy [`templates/Design/`](templates/Design/README.md) and fill in its inputs and provenance before running solvers.

## Repository layout

| Location | Responsibility |
|---|---|
| `pattern design/` | Shared `QComponent`s ([inventory](skill/qcomponents.md)) |
| `code/` | [Shared utilities, circuit models, and run recording](code/README.md) |
| `skill/` | Design-agnostic SOPs and acceptance criteria ([index](skill/skill.md) · [EM defaults](skill/em_simulation.md) · [GDS export](skill/gds_export.md) · [fab deliveries](skill/fabrication_output.md)) |
| `note/` | General physics, formulas, and methods ([Purcell filter](note/purcell_filter_method.md) · [tunable-coupler CZ](note/tunable_coupler_cz_method.md) · [FTF architecture](note/ftf_workflow.md) · [NTHU TiAu alignment](note/tiau_marker_alignment.md)) |
| `Design*/` | Self-contained chip designs (`layout/`, `simulation/`, `results/`, `parameters/`, `gds/`, `fig/`, `notes/`, `design.json`); discovered automatically via `project_paths.all_design_dirs()` |
| `Fabrication output/` | [Dated fabrication deliveries and e-beam marker templates](Fabrication%20output/README.md) |
| `templates/Design/` | [Standard template for new Designs](templates/Design/README.md) |

> **Separation rule:** `skill/` and `note/` hold general, design-agnostic rules and methods only. All chip-specific numbers, solver overrides, GDS links, and tuning history belong inside each `Design*/` folder.

## Workflow and provenance

```text
Targets & constraints → geometry & quantum inputs → Q3D / HFSS → circuit extraction → scqubits / dynamics → confirmed GDS & fab delivery
```

- **Input boundaries:** `EJ` / `EL` are user or quantum-model inputs, and HFSS lumped `L` is an EM input—never label them as direct HFSS outputs.
- **Directory convention:** each Design classifies code and results under `{general, qubit, resonator, flux}` (creating only categories used); new runs record `run.json`, `code/`, `raw/`, and `summaries/` / `figures/` / `logs/` via `code/general/run_record.py` ([organization note](note/project_organization.md)).
- **Version adoption:** select adopted runs explicitly in `design.json` rather than by newest timestamp. After any geometry change, earlier results belong to the earlier revision and must be marked `needs_revalidation`.
- **Release confirmation:** before committing/pushing final parameters or releasing fabrication masks, complete the confirmation checklist in [skill/parameter_provenance.md](skill/parameter_provenance.md).
