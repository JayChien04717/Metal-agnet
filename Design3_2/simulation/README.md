# Design3_2 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [prepare_fx1_port_review.py](flux/prepare_fx1_port_review.py) | flux |
| [run_remaining.py](general/run_remaining.py) | general |
| [run_targeted.py](general/run_targeted.py) | general |
| [run_verification.py](general/run_verification.py) | general |
| [update_readme.py](general/update_readme.py) | general |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
