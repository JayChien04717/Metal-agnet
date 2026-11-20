# DesignNN: enter the design name

> Copy the entire templates/Design directory to the project root as DesignNN, then replace placeholders in this document and JSON files. This template contains no executable layout or simulation code; design code must explicitly load JSON. Null means unset, not zero.

## 1. Design purpose and current status

Describe topology, qubit types, chip size, and targets. Status: draft. Specify geometry_revision and adopted results in [design.json](design.json); explicitly state when simulations have not been run.

## 2. Layout and artifacts

After generation, add full-chip/local figures and GDS links, identifying simulation or fab, geometry revision, and source run. Do not create links to artifacts that do not exist.

## 3. Parameter summary and provenance

| Parameter | Target | Adopted value | Unit | Source/status |
|---|---|---|---|---|
| TBD | TBD | Unset | TBD | draft |

[Geometry inputs](parameters/geometry.json) · [Quantum inputs](parameters/quantum.json) · [Solver settings](parameters/solver.json) · [Provenance](parameters/provenance.json)

EJ/EL come from user settings or explicit quantum-target solutions; fixed values must not be automatically retargeted. Record derivation/fitting methods for EC and coupling. Reference exact values from designated runs rather than maintaining competing official copies.

## 4. Simulation results and applicability

No adopted results yet. For each completed result, list run ID, raw-data/analysis links, geometry revision, validation evidence, and limitations. Update adopted_results in design.json; do not select by file timestamp.

## 5. Reproduction entry points, environments, and Design-specific overrides

No executable entry point yet. After implementing `layout/build.py` and required simulation scripts, record the working directory, complete commands, Python environment, and tool versions. Complete inputs before execution and use a unique ID for every run.

**Design-specific rules and solver overrides:** document any deliberate deviation from the repository defaults in `skill/em_simulation.md` or `skill/gds_export.md` here (for example, reduced HFSS driven pass budgets on large multi-channel chips, custom port assignments, or special landing-pad options) and record the same override in `run.json`.

## 6. Outstanding work and directories

- Fill in design inputs, provenance, and fixed/adjustable constraints.
- Implement layout and required solver workflows; document simulation/fab differences.
- Complete DRC, output-consistency, EM-convergence, and quantum-truncation checks.

```text
DesignNN/
  README.md
  design.json
  parameters/       geometry.json, quantum.json, solver.json, provenance.json
  layout/           Chip assembly; components/ holds design-specific QComponents
  simulation/       general/、qubit/、resonator/、flux/
  results/          Classified by purpose; per-run run.json, code, raw data, and summaries
  gds/              Exported artifacts
  fig/              Structure figures
  notes/            Design tradeoffs, tuning history, and exceptions
```

Shared rules are in root AGENTS.md and skill/parameter_provenance.md. Historical designs may retain paths, but design.json must explicitly record their mapping.

## Confirmation before commit/push

Prepare the final candidate, parameter snapshot, runs, GDS/figures, and unvalidated items for user confirmation; record it in design.json. Candidate selection and complete-parameter confirmation are separate; subsequent adjustments require renewed confirmation. README presents only the last confirmed version. Move old values to notes and leave unvalidated new results pending.
