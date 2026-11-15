# Superconducting Quantum Chip Layout Project

Build layouts with Qiskit Metal, obtain electromagnetic data with Q3D/HFSS, construct quantum models with scqubits, and analyze dynamics using QuTiP or explicitly documented numerical methods. Each Design records the tools actually used.

## Getting started

1. Read the [AI working rules](AGENTS.md) and [parameter provenance rules](skill/parameter_provenance.md).
2. Select a Design and check its current revision, result applicability, and execution environment.
3. For a new design, copy the [Design template](templates/Design/README.md) and complete its inputs and provenance before execution.

## Where content belongs

| Location | Responsibility |
|---|---|
| `pattern design/` | Shared QComponents; GDS tools have moved to code/general/ |
| `code/` | [Shared utilities and run recording](code/README.md) |
| `skill/` | Required workflows and acceptance criteria; [rules index](skill/skill.md) |
| `note/` | General principles, calculation methods, tips, and references |
| `Design*/` | Design inputs, code, simulation results, GDS, figures, and history |
| `templates/Design/` | Standard template for new Designs; unset values do not constitute an executable design |

Skills describe procedures and acceptance criteria; notes explain principles. Keep design values and tradeoffs in the corresponding Design. README is a summary; link exact results to a specific run.

## Design index

[Design1](Design1/README.md) · [Design2](Design2/README.md) · [Design3](Design3/README.md) · [Design3_2](Design3_2/README.md) · [Design3_3](Design3_3/README.md) · [Design4](Design4/README.md) · [Design5](Design5/README.md) · [Design6: organization example](Design6/README.md)

## Workflow and provenance

Design targets/user constraints → geometry and quantum inputs → Q3D/HFSS → circuit extraction/fitting → scqubits → dynamics analysis.

EJ and EL are specified by the user or obtained through quantum-model design with explicit targets and constraints; they are not direct HFSS outputs. EM lumped L is an input. After geometry changes, earlier results still belong to the earlier revision and do not automatically validate the new one.

## Directory structure and compatibility

All Designs use `layout/`, `simulation/{general,qubit,resonator,flux}/`, and `results/{general,qubit,resonator,flux}/`. Create only categories that are used. Mixed historical data remain in general/legacy; consult each Design's results/index.csv.

[Organization record](note/project_organization.md) · [New Design template](templates/Design/README.md). Existing parameters and geometry were not changed. design.json indexes adopted versions; it does not automatically configure every legacy script.

Before commit/push, complete the final version and parameter confirmation required by AGENTS.md. README presents the confirmed version; older data remain in historical records.
