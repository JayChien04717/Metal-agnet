# Design4 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [analyze6.py](resonator/analyze6.py) | resonator |
| [cai_fit_plot.py](resonator/cai_fit_plot.py) | resonator |
| [hfss6.py](resonator/hfss6.py) | resonator |
| [hfss6_meshlink.py](resonator/hfss6_meshlink.py) | resonator |
| [hfss_hammer_driven.py](resonator/hfss_hammer_driven.py) | resonator |
| [hfss_hammer_straps.py](resonator/hfss_hammer_straps.py) | resonator |
| [interp_sweep.py](resonator/interp_sweep.py) | resonator |
| [modal_capacitance.py](resonator/modal_capacitance.py) | resonator |
| [q3d6.py](resonator/q3d6.py) | resonator |
| [q3d_hammer_gap.py](resonator/q3d_hammer_gap.py) | resonator |
| [run_linked_gui.py](resonator/run_linked_gui.py) | resonator |
| [s21_model_vs_eigen.py](resonator/s21_model_vs_eigen.py) | resonator |
| [verify_model.py](resonator/verify_model.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
