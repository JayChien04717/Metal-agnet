# Design3 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [hfss_mutual_fx.py](flux/hfss_mutual_fx.py) | flux |
| [plot_mutual_z.py](flux/plot_mutual_z.py) | flux |
| [run_mutual_current.py](flux/run_mutual_current.py) | flux |
| [run_mutual_d3.py](flux/run_mutual_d3.py) | flux |
| [sweep_mutual.py](flux/sweep_mutual.py) | flux |
| [q3d_ftf.py](qubit/q3d_ftf.py) | qubit |
| [sweep_FX.py](qubit/sweep_FX.py) | qubit |
| [sweep_FX_charge_line.py](qubit/sweep_FX_charge_line.py) | qubit |
| [sweep_TC.py](qubit/sweep_TC.py) | qubit |
| [zz_analysis.py](qubit/zz_analysis.py) | qubit |
| [hfss_eigen_res.py](resonator/hfss_eigen_res.py) | resonator |
| [q3d_kappa.py](resonator/q3d_kappa.py) | resonator |
| [q3d_tee.py](resonator/q3d_tee.py) | resonator |
| [scq_dispersive.py](resonator/scq_dispersive.py) | resonator |
| [sweep_readout_pad.py](resonator/sweep_readout_pad.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
