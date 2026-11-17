# Design6: tunable Xmon–coupler–Xmon chip

**Current wafer layout 2026-11-17_wafer_r10.** This revision uses the four tail lengths from the trimmed loaded HFSS run. The readout CPW pairs retain 200 µm centre spacing and 166 µm etched inner-edge spacing. The die includes four T test junctions and the ±2.00 mm side launches. Full-chip electrical validation is pending.

[Fab GDS](gds/Design6_2026-11-17_wafer_r10_fab.gds) · [Simulation geometry GDS](gds/Design6_2026-11-17_wafer_r10_sim.gds) · [Fab mask preview](fig/Design6_2026-11-17_wafer_r10_fab_check.png) · [r9 parent mask](gds/Design6_2026-11-16_wafer_r9_fab.gds) · [Wafer index](../note/wafer_placement_2026-11-16.md)

## Loaded readout S21

![Design6 r10 loaded HFSS S21 across Q2 and Q1 with an orange dashed four-pole fit](fig/Design6_2026-11-17_wafer_r10_loaded_s21.png)

The blue trace covers 6.55–7.25 GHz and includes both Xmons and readout claws. The orange dashed four-pole fit spans the 6.751–7.056 GHz fine sweep. The GDS uses the lengths from this HFSS run.

| Channel | Filter / resonator length (mm) | Bare f_p / f_r (GHz) | f_r − f_p (MHz) | J (MHz) | κ_p/2π (MHz) |
|---|---:|---:|---:|---:|---:|
| Q1 | 3.880590 / 4.053848 | 6.99175 / 6.98892 | −2.83 | 15.96 | 36.15 |
| Q2 | 3.993377 / 4.137843 | 6.80975 / 6.82477 | +15.02 | 15.25 | 31.66 |

The HFSS fixture reached ΔS = 0.070565 after three adaptive passes, above its 0.01 stop target. The two hybrid linewidths are 16.15/20.00 MHz (Q1) and 23.48/8.17 MHz (Q2), in ascending frequency order. The GDS export used the existing simulation evidence. [Raw S21](results/resonator/20261117_193641_readout_loaded_r9_trim_98139f87/raw/S_fine.csv) · [Fit table](results/resonator/20261117_201559_compare_loaded_trim_r9_5460120e/fit_comparison.csv) · [Before-and-after comparison](notes/readout_s21_20261117.md#loaded-hfss-tail-trim-candidate)

The etched die ID `TXT00` is 50 µm inside the top and left die edges. The 5 × 5 µm ground perforations use a 20 µm centre pitch and a 547.5 µm four-edge keep-out. See the [wafer cheese rule](../skill/gds_export.md#3-checks-after-an-export).


> **Directory reorganization (2026-11-15)**: [Layout code](layout/) · [Simulation entry points](simulation/README.md) · [Results index](results/README.md) · [Adopted version](design.json). Solvers were not rerun. Moving historical values does not make them final confirmed values. Create a new run for new executions and confirm final parameters before commit/push.


## 1. Design purpose and current status

A 10 × 10 mm chip with two asymmetric-SQUID Xmon qubits (Q1, Q2), a central tunable SQUID transmon coupler (TC), and Purcell-filtered readout for each qubit.

The wafer revision retains the **B2 flux-line geometry** selected on 2026-11-15. The complete parameter snapshot awaits confirmation before commit/push. The 2026-11-14 core, readout and CZ records are historical revision-A results.

[Design index](design.json) · [Parameter provenance](parameters/provenance.json) · [Full historical report](notes/report_2026-11-15.md) · [Parameter overview (all current values with status, 2026-11-17)](notes/parameters_overview_20261117.md)

## 2. Earlier B2 layout and artifacts

![Full chip](fig/Design6_2026-11-15_layout.png)
![Core structure](fig/Design6_2026-11-15_core.png)

- [Simulation GDS](gds/previous/Design6_2026-11-15_sim.gds): simulation electrodes, with junction markers on layer 50.
- [Fab GDS](gds/previous/Design6_2026-11-15_fab.gds): fabrication contacts, fillets, and ground perforation; distinguish this from simulation geometry.
- [Assembly entry point](layout/design6_layout.py); dedicated components are xmon.py and claw_coupler.py in the same directory.

The images in this section document the November 15 B2 layout. The current wafer files and mask preview are linked above.

## 3. Final candidate B2 parameters

| Item | B2 setting | Source/status |
|---|---|---|
| Flux line | L shape, sep 9 µm, bar 45 µm, fbl_x 12 µm | Existing B2 report; candidate confirmed by the user |
| EJ/EL, EC, coupling | Complete final parameter set awaits confirmation | Do not reuse revision-A results as official B2 values |

Confirming the B2 candidate does not automatically adopt historical retargeted EJ values. EJ/EL must still follow user constraints or quantum-model settings; HFSS does not directly provide EJ/EL. Historical values are in the [earlier report](notes/report_2026-11-15.md) and [historical provenance index](parameters/provenance.json).

## 4. B2 result status

- Local flux: the existing report records B2 results; raw-run mapping remains incomplete, so these are not included in the official numerical table.
- Launch pads: local [pad_s11 data](results/general/2026-11-15/pad_s11) exist; complete parameters and correspondence to B2 require final verification.
- Core Q3D on r9: the user adopted pass 20 as the current design input. AEDT reported Δ = 0.27546% against its 0.1% target. [r9 and historical v8 comparison](notes/parameters_overview_20261117.md#core-q3d-on-r9)
- Readout and XY: [r10 loaded readout](#loaded-readout-s21) · [earlier readout-only and Q3D model comparisons](notes/readout_s21_20261117.md) · [local XY fixture](notes/parameters_overview_20261117.md#xy-termination-candidates).
- CZ: [r9 coherent comparison using the adopted Q3D core](notes/cz_coupler_flux_20261117.md#r9-hamiltonian). The r9 HFSS core fit and full-device validation remain pending.

## 5. Reproduction entry points and environments

Start from the project root. These existing script entry points do not guarantee reproduction of historical revision-A geometry; establish version correspondence when code snapshots are missing. Use a new date/tag and check for existing output before execution.

Wafer GDS reproduction uses the current `wafer_candidate` in [design.json](design.json) and reads back the selected regions with KLayout. Choose unused output filenames.

```powershell
conda run -n qmetal python -B "code/general/export_wafer_test_junctions.py" --design Design6 --output Design6_r10_review_fab.gds
conda run -n qmetal python -B "code/general/export_wafer_test_junctions.py" --design Design6 --reference-only --output Design6_r10_review_sim.gds
```

Core Q3D for r9 (qmetal, AEDT 2021 R2):

```powershell
conda run -n qmetal python -X utf8 "Design6/simulation/qubit/q3d_core_r9.py" --dry-run
conda run -n qmetal python -X utf8 "Design6/simulation/qubit/q3d_core_r9.py"
```

Add B2 quantum-analysis commands after final EJ/EL constraints, parameter files, and corresponding EM runs are established. See the [earlier report §8](notes/report_2026-11-15.md) for historical commands.

## 6. Outstanding work and directory compatibility

- Complete B2 flux raw-run mapping and execution snapshots.
- Check the effects of launch-pad/B2 geometry on core, readout, and CZ results.
- Record drawn Lj, quantum-retargeted Lj, and fabrication values separately.
- Mark results validated only after recording convergence and truncation evidence.
- Coupler-flux-only CZ exploration (2026-11-17, scqubits, hypothetical targets, not adopted): [note](notes/cz_coupler_flux_20261117.md). The current targets have no ZZ-off point (Δ12 > |α|); a straddling set (Δ12 ≈ 150 MHz, g_qc ≈ 61 MHz, g12 ≈ 4 MHz, coupler f_max 6.8 GHz) gives idle ZZ +5 kHz and a 160 ns CZ with coherent F 0.9995.

This update reorganized code and result paths without rerunning solvers. Layout entry points are in layout/, quantum analysis in simulation/qubit/, and results are classified as general/qubit/resonator/flux. design.json and provenance.json remain provenance/adoption indexes; they do not yet control all parameter loading in legacy scripts.

## Appendix A. Input capacitor and readout margin

The series input capacitor \(C_\mathrm{in}\) reflects part of the incident signal at P7 and redirects readout emission toward P8. The model uses a lossless capacitor and matched 49.9 Ω feedline arms. The input insertion loss is evaluated at 6.9 GHz and comes from reflection. The P8 share is the minimum across the four r10 hybrid-mode frequencies, with the capacitor position optimized separately for each value. [Heinsoo et al., Appendix C](https://arxiv.org/pdf/1801.07904#page=8)

| \(C_\mathrm{in}\) (fF) | Input insertion loss at 6.9 GHz (dB) | Minimum share of feedline emission at P8 |
|---:|---:|---:|
| 23 | 20.1 | 99.1% |
| **73** | **10.4** | **92.2%** |
| 140 | 5.7 | 78.5% |

The preferred 73 fF candidate exceeds the 75% output target in the model at an electrical distance of about 2.42 mm from the Q2 tee toward P7. The estimate omits capacitor parasitics and changes to the mode frequencies and linewidths. This candidate remains pending layout and EM validation.

In the r10 loaded HFSS baseline, the bare filter and resonator frequencies differ by 2.83 MHz for Q1 and 15.02 MHz for Q2. The Q2 detuning is comparable to its fitted \(J\) of 15.25 MHz. A scqubits calculation using the adopted r9 core Q3D parameters and fixed layout \(L_J\), combined with the r10 fitted poles, estimates the full \(|0\rangle\)-to-\(|1\rangle\) frequency separation of the hybrid modes. [Core Q3D inputs](results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/summary_out.json) · [Dispersive calculation](simulation/qubit/readout_qubit_r9.py) · [Loaded fit summary](results/resonator/20261117_201559_compare_loaded_trim_r9_5460120e/summary.json)

| Channel | Hybrid mode with larger shift | Full state separation (MHz) | Fitted FWHM (MHz) | Separation / FWHM |
|---|---|---:|---:|---:|
| Q1 | Lower frequency | 1.22 | 16.15 | 0.075 |
| Q2 | Higher frequency | 1.20 | 8.17 | 0.147 |

The estimated shifts support further readout validation. Q1 has less spectral contrast because its preferred mode is broader. Assessing 99% single-shot assignment fidelity requires the state-dependent output IQ separation at the intended readout power and integration time, along with the measurement-chain efficiency and qubit relaxation. After adding \(C_\mathrm{in}\), the reflected loading can shift \(f_p\) and change both hybrid linewidths. The corresponding Purcell-limited lifetimes still require decay-rate estimates at the qubit frequencies of 4.81 GHz and 4.41 GHz. No solver was rerun for this appendix. [Hybrid-mode readout analysis](https://doi.org/10.1103/PRXQuantum.5.040326) · [Readout-system characterization](https://arxiv.org/abs/2402.00413)
