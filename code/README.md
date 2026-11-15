# Shared code

| Category | Contents |
|---|---|
| general | GDS, plotting, AEDT connections, path compatibility, and run recording |
| qubit | Shared capacitance-matrix and quantum-analysis algorithms |

Create resonator/flux categories when genuinely shared code is available. Scripts tied to Design-specific constants remain within their Design. Shared QComponents stay in the root pattern design/ directory.

`project_paths.py` and `path_migration.json` resolve old data paths without duplicating old directories or creating shortcuts. Historical manifests remain unchanged; use the index to find their current locations.

## New runs

Run from the project root (prepares a record only; does not start a solver):

```powershell
conda run -n qmetal python code/general/run_record.py Design3_3 flux B2_check --entrypoint Design3_3/simulation/flux/recheck_readout_flux.py
```

The command returns a unique run directory. Pass confirmed inputs with --input and parent results with --upstream. The solver entry point must explicitly accept this directory; creating a record alone does not redirect output away from old results.

In Python, call `create_run(..., parameters=actual_settings)`, write outputs into the returned directory, then call `finish_run(..., status='completed', summary=results)`; record failures as failed. Keep validation unknown until convergence is checked. Record actual geometry_revision and solver_version values rather than guessing.

To cover cross-Design dependencies, the code snapshot includes all active Python files and the path index, excluding copies inside historical results. Third-party package versions are recorded; this does not claim to archive the complete execution environment.
