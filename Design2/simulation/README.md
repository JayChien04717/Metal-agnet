# Design2 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [hfss_mutual_fx.py](flux/hfss_mutual_fx.py) | flux |
| [run_mutual.py](flux/run_mutual.py) | flux |
| [analyse.py](general/analyse.py) | general |
| [q3d_ftf.py](qubit/q3d_ftf.py) | qubit |
| [run_q3d.py](qubit/run_q3d.py) | qubit |
| [transition_estimates.py](qubit/transition_estimates.py) | qubit |
| [hfss_eigen_res.py](resonator/hfss_eigen_res.py) | resonator |
| [linewidth.py](resonator/linewidth.py) | resonator |
| [modal_linewidth.py](resonator/modal_linewidth.py) | resonator |
| [plot_readout_modes.py](resonator/plot_readout_modes.py) | resonator |
| [q3d_kappa.py](resonator/q3d_kappa.py) | resonator |
| [readout_local_eigen.py](resonator/readout_local_eigen.py) | resonator |
| [run_linear_readout.py](resonator/run_linear_readout.py) | resonator |
| [run_readout.py](resonator/run_readout.py) | resonator |
| [run_res_modes.py](resonator/run_res_modes.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
