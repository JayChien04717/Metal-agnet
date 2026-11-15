# Simulation results

Create general/qubit/resonator/flux categories only when results exist. Each run retains run.json, code/, raw/, and any needed summary.json, figures, and logs. Use shared code/general/run_record.py for a unique directory; entry points must explicitly write there.

Reanalysis creates a new run with upstream references; never overwrite raw data. Select adopted results explicitly in design.json.
