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
| [hfss_core.py](qubit/hfss_core.py) | qubit |
| [hfss_pad.py](qubit/hfss_pad.py) | qubit |
| [plot_ej_layout.py](qubit/plot_ej_layout.py) | qubit |
| [q3d_core.py](qubit/q3d_core.py) | qubit |
| [q3d_core_r9.py](qubit/q3d_core_r9.py) | qubit |
| [epr_readout.py](resonator/epr_readout.py) | resonator |
| [hfss_readout.py](resonator/hfss_readout.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
