# Design5_6_purcell-res simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [modes5.py](resonator/modes5.py) | resonator |
| [purcell_resonator_pair_circuit_model.py](resonator/purcell_resonator_pair_circuit_model.py) | resonator |
| [purcell_resonator_pair_compare_plot.py](resonator/purcell_resonator_pair_compare_plot.py) | resonator |
| [purcell_resonator_pair_s21_reanalyze.py](resonator/purcell_resonator_pair_s21_reanalyze.py) | resonator |
| [runstore.py](resonator/runstore.py) | resonator |

**Removed 2026-11-25** (depended on the deleted `Design4_Single_purcell-res` folder or on the removed `design5_layout` / `chip6_*` layouts; in git history): `analyze5.py`, `chip6_d4_full.py`, `chip6_kappa_env.py`, `chip6_theory.py`, `eig_chip6.py`, `hfss_chan5.py`, `hfss_chip6.py`, `purcell_resonator_pair_validate.py`, `q3d_full5.py`, `q3d_hammer_check.py`, `q3d_hammer_final.py`, `solve_chip6.py`, `theory5.py`.

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
