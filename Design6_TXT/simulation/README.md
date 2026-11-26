# Design6_TXT simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [compare_flux_variants.py](flux/compare_flux_variants.py) | flux |
| [flux_metrics.py](flux/flux_metrics.py) | flux |
| [flux_pulse.py](flux/flux_pulse.py) | flux |
| [flux_resweep.py](flux/flux_resweep.py) | flux |
| [flux_square_spectro.py](flux/flux_square_spectro.py) | flux |
| [flux_tail_scan.py](flux/flux_tail_scan.py) | flux |
| [flux_xtalk_matrix.py](flux/flux_xtalk_matrix.py) | flux |
| [hfss_flux.py](flux/hfss_flux.py) | flux |
| [cz_scqubits.py](qubit/cz_scqubits.py) | qubit |
| [cz_qutip_r9.py](qubit/cz_qutip_r9.py) | qubit (QuTiP check of the r9 CZ: sesolve, truncation, sensitivity, mesolve T1/T2) |
| [hfss_core.py](qubit/hfss_core.py) | qubit |
| [hfss_pad.py](qubit/hfss_pad.py) | qubit |
| [plot_ej_layout.py](qubit/plot_ej_layout.py) | qubit |
| [q3d_core.py](qubit/q3d_core.py) | qubit |
| [q3d_core_r9.py](qubit/q3d_core_r9.py) | qubit |
| [compare_loaded_s21_r9.py](resonator/compare_loaded_s21_r9.py) | resonator |
| [compare_loaded_trim_r9.py](resonator/compare_loaded_trim_r9.py) | resonator |
| [compare_readout_cin.py](resonator/compare_readout_cin.py) | resonator |
| [plot_loaded_s21.py](resonator/plot_loaded_s21.py) | resonator |
| [purcell_resonator_pair_validate.py](resonator/purcell_resonator_pair_validate.py) | resonator |
| [readout_network.py](resonator/readout_network.py) | resonator |
| [readout_poles.py](resonator/readout_poles.py) | resonator |
| [readout_q3d_check.py](resonator/readout_q3d_check.py) | resonator |
| [wallraff_check.py](resonator/wallraff_check.py) | resonator |

**Removed 2026-11-25** (read data from the deleted `Design4_Single_purcell-res` folder or imported a script that did; in git history): `resonator/epr_readout.py`, `hfss_readout.py`, `parallel_hammer_check.py`, `readout_cin_r10.py`, `readout_compare.py`, `readout_length_r9.py`, `readout_loaded_r9.py`, `readout_model.py`, `readout_model_q3d.py`, `readout_q3d_full.py`, `readout_s21.py`. Their saved runs under `results/` remain as historical records.

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
