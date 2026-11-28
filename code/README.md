# Shared code (`code/`)

Shared, design-agnostic utilities organized by function. `code/general/project_paths.py` (`configure(__file__)`) automatically adds all active `code/` subfolders to `sys.path` so bare imports (`from gds_tools import ...`, `from check_keepouts import ...`, `from wafer_cheese_grid import ...`) work uniformly across Designs.

| Folder | Contents |
|---|---|
| `general/` | Shared bootstrap, tracking & PDF reports: `project_paths.py`, `path_migration.json`, `run_record.py`, `sync_tracked_results.py`, `parameter_pdf.py` |
| `aedt/` | Ansys AEDT / HFSS session & file helpers: `aedt_session.py`, `aedt_messages.py`, `hfss_dirs.py`, `prune_aedt_designs.py`, `save_aedt.py` |
| `layout/` | Layout DRC, GDS export & figures: `gds_tools.py`, `check_keepouts.py`, `check_metal_width.py`, `route_length_solver.py`, `plot_design_figs.py` (+ `design_figures.json`) |
| `fab/` | 4-inch wafer assembly, tone inversion, cheese & alignment markers: `assemble_wafer.py`, `assemble_wafer_mask2_grid.py`, `invert_wafer_tone.py`, `wafer_cheese_grid.py`, `wafer_cheese_keepout.py`, `wafer_region_overlay.py`, `marker_ground_opening.py`, `tiau_marker_overlay.py`, `tsri_marker_overlay.py`, `check_marker_references.py` |
| `ebeam/` | E-beam junction placement, checks & plots (in development): `ebeam_junction_placement.py`, `check_ebeam_junctions.py`, `plot_ebeam_junction_parameters.py` |
| `qubit/` | Qubit physics & Q3D analysis: `cap_matrix.py`, `jj_resistance_model.py`, `preview_fluxonium_delft.py` |
| `resonator/` | Purcell / readout resonator models & sim: `cpw_pair_cap2d.py`, `purcell_readout_model.py`, `purcell_resonator_pair_sim.py` |

Scripts tied to Design-specific constants remain inside each `Design<N>/`. Shared `QComponent`s stay in root `pattern design/`.

`aedt_session.py` closes an AEDT that a run launched when the run ends (saves projects, quits, watchdog kill after a hard kill); it is installed by importing `project_paths` and leaves an already-open AEDT alone. Set `AEDT_KEEP_OPEN=1` to keep AEDT open. See `skill/em_simulation.md` §4.

## New runs

Run from the project root (prepares a record only; does not start a solver):

```powershell
conda run -n qmetal python code/general/run_record.py <Design> <category> <label> --entrypoint <script>
```

In Python, call `create_run(..., parameters=actual_settings)`, write outputs into the returned directory, then call `finish_run(..., status='completed', summary=results)`; record failures as failed. Keep validation unknown until convergence is checked.
