# Design3 simulation code

Classify by research purpose rather than solver. general covers cross-purpose checks; HFSS/Q3D/scqubits names remain in entry-point filenames.

| Entry point | Category |
|---|---|
| [recheck_readout_flux.py](flux/recheck_readout_flux.py) | flux |
| [report_readout_flux.py](flux/report_readout_flux.py) | flux |
| [verify.py](general/verify.py) | general |
| [compare.py](qubit/compare.py) | qubit |
| [optimize_pulse.py](qubit/optimize_pulse.py) | qubit |
| [q3d.py](qubit/q3d.py) | qubit |
| [recheck_pulse.py](qubit/recheck_pulse.py) | qubit |
| [run_q3d.py](qubit/run_q3d.py) | qubit |
| [sequence.py](qubit/sequence.py) | qubit |
| [sweep.py](qubit/sweep.py) | qubit |
| [validate_gate.py](qubit/validate_gate.py) | qubit |
| [probe_readout.py](resonator/probe_readout.py) | resonator |

Legacy input/output interfaces are retained; some scripts still write fixed filenames, so specify a new output location before execution. `code/general/run_record.py` creates a run and snapshots code/parameters, but does not silently redirect legacy output. Reclassification changes paths only; it does not imply rerunning or validating results.
