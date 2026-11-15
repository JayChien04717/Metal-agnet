# Design5 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [analyze5.py](resonator/analyze5.py) | resonator |
| [chip6_d4_full.py](resonator/chip6_d4_full.py) | resonator |
| [chip6_kappa_env.py](resonator/chip6_kappa_env.py) | resonator |
| [chip6_theory.py](resonator/chip6_theory.py) | resonator |
| [eig_chip6.py](resonator/eig_chip6.py) | resonator |
| [hfss_chan5.py](resonator/hfss_chan5.py) | resonator |
| [hfss_chip6.py](resonator/hfss_chip6.py) | resonator |
| [modes5.py](resonator/modes5.py) | resonator |
| [q3d_full5.py](resonator/q3d_full5.py) | resonator |
| [runstore.py](resonator/runstore.py) | resonator |
| [solve_chip6.py](resonator/solve_chip6.py) | resonator |
| [theory5.py](resonator/theory5.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
