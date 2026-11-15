# Design6: tunable Xmon–coupler–Xmon chip


> **Directory reorganization (2026-11-15)**: [Layout code](layout/) · [Simulation entry points](simulation/README.md) · [Results index](results/README.md) · [Adopted version](design.json). Solvers were not rerun. Moving historical values does not make them final confirmed values. Create a new run for new executions and confirm final parameters before commit/push.


## 1. Design purpose and current status

A 10 × 10 mm chip with two asymmetric-SQUID Xmon qubits (Q1, Q2), a central tunable SQUID transmon coupler (TC), and Purcell-filtered readout for each qubit.

The user selected **B2 as the final candidate** on 2026-11-15. The layout is **2026-11-15 B2 flux lines + updated launch pads**. The B2 candidate is confirmed; the complete parameter snapshot awaits confirmation before commit/push. **Only some modules have corresponding checks; the 2026-11-14 core/readout/CZ results remain historical revision-A model results and do not validate the complete B2 chip.**

[Design index](design.json) · [Parameter provenance](parameters/provenance.json) · [Full historical report](notes/report_2026-11-15.md)

## 2. Layout and artifacts

![Full chip](fig/Design6_2026-11-15_layout.png)
![Core structure](fig/Design6_2026-11-15_core.png)

- [Simulation GDS](gds/Design6_2026-11-15_sim.gds): simulation electrodes, with junction markers on layer 50.
- [Fab GDS](gds/Design6_2026-11-15_fab.gds): fabrication contacts, fillets, and ground perforation; distinguish this from simulation geometry.
- [Assembly entry point](layout/design6_layout.py); dedicated components are xmon.py and claw_coupler.py in the same directory.

Figures and GDS are indexed by existing filenames and were not regenerated in this update. Source hashes from the original exports remain unverified.

## 3. Final candidate B2 parameters

| Item | B2 setting | Source/status |
|---|---|---|
| Flux line | L shape, sep 9 µm, bar 45 µm, fbl_x 12 µm | Existing B2 report; candidate confirmed by the user |
| EJ/EL, EC, coupling | Complete final parameter set awaits confirmation | Do not reuse revision-A results as official B2 values |

Confirming the B2 candidate does not automatically adopt historical retargeted EJ values. EJ/EL must still follow user constraints or quantum-model settings; HFSS does not directly provide EJ/EL. Historical values are in the [earlier report](notes/report_2026-11-15.md) and [historical provenance index](parameters/provenance.json).

## 4. B2 result status

- Local flux: the existing report records B2 results; raw-run mapping remains incomplete, so these are not included in the official numerical table.
- Launch pads: local [pad_s11 data](results/general/2026-11-15/pad_s11) exist; complete parameters and correspondence to B2 require final verification.
- Core/readout/CZ: saved revision-A results are historical references only. Official B2 results await validation; do not fill gaps with old values.

## 5. Reproduction entry points and environments

Start from the project root. These existing script entry points do not guarantee reproduction of historical revision-A geometry; establish version correspondence when code snapshots are missing. Use a new date/tag and check for existing output before execution.

Layout export (qmetal; updates the specified date's GDS/figures; GUI requires a desktop; replace YYYY-MM-DD):

```powershell
conda run -n qmetal python "Design6/layout/design6_layout.py" --date YYYY-MM-DD --test-junctions no
```

Add B2 quantum-analysis commands after final EJ/EL constraints, parameter files, and corresponding EM runs are established. Historical commands remain in the [earlier report §8](notes/report_2026-11-15.md), not as official B2 reproduction instructions.

## 6. Outstanding work and directory compatibility

- Complete B2 flux raw-run mapping and execution snapshots.
- Check the effects of launch-pad/B2 geometry on core, readout, and CZ results.
- Record drawn Lj, quantum-retargeted Lj, and fabrication values separately.
- Mark results validated only after recording convergence and truncation evidence.

This update reorganized code and result paths without rerunning solvers. Layout entry points are in layout/, quantum analysis in simulation/qubit/, and results are classified as general/qubit/resonator/flux. design.json and provenance.json remain provenance/adoption indexes; they do not yet control all parameter loading in legacy scripts.
