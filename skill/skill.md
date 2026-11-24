# Qiskit Metal Coding & Layout Workflow

Start here: [project README](../README.md), [AI rules](../AGENTS.md), [parameter provenance](parameter_provenance.md), [new Design template](../templates/Design/README.md).

This document defines the coding conventions for custom Qiskit Metal planar `QComponent`s (`pattern design/`) and the chip layout script workflow (`Design<N>/layout/`).

## Where things live

`skill/` and `note/` hold **general, design-agnostic** rules and methods only. Everything specific to one Design (numbers, tuning logs, run records, Design-specific solver overrides or exceptions) lives in that Design's folder.

| Folder / file | Content |
|---|---|
| `skill/skill.md` | This file: `QComponent` coding rules and layout-script workflow |
| `skill/layout_structure.md` | Three-file `Design<N>/layout/` structure (`placement.json`, `routes.json`, `build.py`), keep-outs, length solver, and fab export flow |
| `skill/em_simulation.md` | Python/AEDT environments, Q3D/HFSS pass & mesh defaults, wire bonds, HFSS lessons |
| `skill/parameter_provenance.md` | Parameter roles, source boundaries, per-run records, adoption, final confirmation |
| `skill/qcomponents.md` | Shared `QComponent` inventory (`pattern design/`), import paths, 4 µm metal-width rule |
| `skill/chip_drc.md` | Layout DRC before GDS export |
| `skill/gds_export.md` | Single-die GDS export: `jj_metal=False`, fab hooks, launch pad standard, cheesing, test coupons, checks |
| `skill/fabrication_output.md` | 4-inch wafer assembly & fab deliveries: `Fabrication output/<YYYYMMDD>/`, streets, die IDs, positive + etch, e-beam markers |
| `skill/ebeam_junction.md` | E-beam junction masks: EJ → Rn → area → mask width, `FluxoniumJunction` / `TransmonJunction`, placement and checks |
| `skill/layout_figures.md` | Standard layout figures rendered from exported GDS |
| `skill/cap_matrix_analysis.md` | Procedure: Q3D capacitance matrix → Ec, J, g, β, κ |
| `skill/mutual_inductance_analysis.md` | Procedure: flux line → loop M from HFSS |
| `skill/purcell_resonator_design.md` | Purcell filter + readout resonator channel (`PurcellReadout`) |
| `skill/two_qubit_gate_design.md` | Two-qubit flux CZ gate design (transmon–TC–transmon) |
| `note/` | General physics and methods (`purcell_filter_method.md`, `tunable_coupler_cz_method.md`, `ftf_workflow.md`, `tiau_marker_alignment.md`, `project_organization.md`) |
| `code/` | Shared utilities and models (`code/qubit/cap_matrix.py`, `code/layout/gds_tools.py`, `code/resonator/purcell_readout_model.py`, `code/general/run_record.py`) |
| `pattern design/` | Shared `QComponent`s (`fluxonium/`, `transmon/`, `resonator/`, `ebeam_junction/`, `marker/`) |
| `Design<N>/` | Design-specific layout, parameters, solver overrides, simulation runs, GDS, figures, and `notes/` |

---

## 1. Custom `QComponent` coding rules

1. **Class structure:** inherit from `QComponent` (`qiskit_metal.qlibrary.core`). Define `default_options = Dict(...)` (including `pos_x`, `pos_y`, `orientation`, `layer`), `component_metadata = Dict(short_name=..., _qgeometry_table_poly=..., _qgeometry_table_path=..., _qgeometry_table_junction=...)`, and `make(self)`.
2. **Parameter parsing (`self.p`):** inside `make()`, always read `p = self.p` (floats in meters) or `self.parse_value(...)`—never do arithmetic on raw string values in `self.options`.
3. **Transform order (rotate then translate):** draw all Shapely geometries (`draw.rectangle`, `draw.Polygon`, `draw.LineString`, rounded corners via `.buffer(-r, join_style=1).buffer(r, join_style=1)`) in local coordinates centered at `(0, 0)`, then apply mirror (if any) → `draw.rotate(geom, p.orientation, origin=(0, 0))` → `draw.translate(geom, p.pos_x, p.pos_y)`.
4. **QGeometry registration (`self.add_qgeometry`):**
   - `'poly'`: solid metal pads/islands (`subtract=False`, default) and etched ground-plane cutouts (`subtract=True`, pocket larger than pad by at least `connection_gap`).
   - `'path'`: CPW traces (`LineString` / `MultiLineString`, requires `width`).
   - `'junction'`: Josephson junctions / SQUIDs (`LineString`, requires `width`); both endpoints must touch/overlap solid metal pads so EM lumped ports excite properly. Never register junctions in `'path'`.
5. **Clean boolean results:** `draw.subtract` can leave zero-area slivers or self-intersecting polygons; drop slivers (`_clean`) before `add_qgeometry` so Ansys does not fail with *Cannot create rectangle with zero length/width*.
6. **Pins (`self.add_pin`):** define a perpendicular trace-width `LineString` and a non-zero **outward** unit normal vector. Rotate and translate the `LineString`; rotate only the normal vector (and swap point order when mirrored so the normal stays outward). Pass `input_as_norm=True` with `width` and `gap`.
7. **Junction orientation & simulation models:**
   - Draw all junction-table lines of a device in the same direction (Dolan bridges require parallel junctions).
   - Provide `jj_metal` (default `'True'`) on any component with junction-related metal so `jj_metal='False'` yields a junction-free mask ([gds_export.md](gds_export.md)).
   - Use **one simulation model per qubit type** across Q3D, HFSS eigenmode, and mutual inductance, matching the mask's pads/stems/pocket/lines. For fluxonium, use `FluxoniumArrayJJ` (`pattern design/fluxonium/fluxonium_nthu.py`: short centre JJ gap + U-shaped two-side lumped-L array joined by a metal closing end, Kron-reduced in Q3D by `code/qubit/cap_matrix.py`)—never a single straight inductor.
   - Draw shorted flux lines as ground for Q3D (`flux_as_ground='True'`) to avoid net-overlap errors.

### `QComponent` review checklist

- [ ] All arithmetic in `make()` uses `self.p` (or `self.parse_value`), never raw `self.options`.
- [ ] Shapes are drawn locally, then rotated before translating; zero-area slivers cleaned.
- [ ] Ground cutouts use `subtract=True`; solid metal uses default `subtract=False`.
- [ ] Every pin normal points outward (`design.components[name].pins`) and rotates/mirrors properly.
- [ ] Junctions are in the `'junction'` table, touch metal pads at both ends, and are parallel.
- [ ] Ground strips satisfy the 4 µm minimum metal-width rule (`metal_width.require_strip()`, [qcomponents.md](qcomponents.md)).

---

## 2. Layout script workflow & stage gates

Each Design builds its chip under `Design<N>/layout/` following [layout_structure.md](layout_structure.md) (`placement.json`, `routes.json`, and `build.py` exposing `build_design(variant='sim'|'fab')`; no notebook). GDS export uses `code/layout/gds_tools.py` ([gds_export.md](gds_export.md)); standard figures are rendered from the clean reference GDS (`*_sim.gds`) by `code/layout/plot_design_figs.py` ([layout_figures.md](layout_figures.md)).

### 2.0 Three-stage gate rule (do not jump ahead)
1. **Stage 1 — Draft layout & visual review:**
   - Implement only the simulation `QComponent`(s) and `layout/` (`placement.json`, `routes.json`, `build.py`), placing the target e-beam site's alignment marker component (e.g. `TSRIA12Mark` or `TiAuChipMarker`, when applicable) directly in `placement.json` so marker openings appear in `fig/chip.png`.
   - Verify basic layout DRC: no short-segment warnings (`≥ 2 × R`), all pin `middle` and outward `normal` vectors, filleted resonator lengths equal to targets, no crossing routes, keep-out clearance (`check_keepouts`, ≥ 100 µm from **both** NTHU TiAu and TSRI A12 marker openings to active routes/circuitry/coupons), and component-level 4 µm minimum metal width (`metal_width.require_strip()` / `check_component`).
   - Export **only** the clean reference `*_sim.gds` needed to render `fig/chip.png` and `fig/<part>.png` via `code/layout/plot_design_figs.py`.
   - **Gate 1 (STOP):** Show the layout figures and basic DRC summary to the user and wait for feedback. **Do not** build `*_fab.py` Dolan hooks, place `FluxoniumJunction` / `TransmonJunction` e-beam layers, add corner test coupons or cheesing, export `*_fab.gds`, or run mask-level fab output checks unless the user explicitly asks for them.
2. **Stage 2 — EM simulation & parameter tuning:**
   - Only after the user approves the draft layout (or explicitly asks to run simulations), run Q3D / HFSS ([em_simulation.md](em_simulation.md)) and iterate geometry until electrical targets are confirmed.
3. **Stage 3 — Fab mask, e-beam junctions & wafer export:**
   - Only after simulation parameters are confirmed and the user requests fabrication output, follow [gds_export.md §2.1](gds_export.md#21-fab-gds-workflow) (ask the checklist questions, add `*_fab` / e-beam layers, export `*_fab.gds`, and run full mask checks).

### 2.1 Knobs, paths, and derivation
* **Never hardcode the Design folder name** inside layout or simulation scripts. Resolve the enclosing Design folder dynamically via `DESIGN = project_paths.design_dir(__file__)` (and pass `DESIGN.name` to `run_record.create_run(...)`), so renaming a Design folder never breaks its scripts. Shared tools under `code/` should discover Designs via `project_paths.all_design_dirs()` (directories containing `design.json`).
* Keep all tunable dimensions in module-level dicts/constants at the top of the script and accept overrides in `build_design(...)` so sweeps never edit source code.
* Derive dependent geometry (e.g., coupler dimensions from qubit pad dimensions) instead of duplicating constants.

### 2.2 Routing and CPW spacing
* Use `RouteAnchors` with explicit anchors for control lines (Pathfinder fails easily on crowded chips), `RouteStraight` between neighbouring tees, and deterministic serpentines for long resonators.
* `RouteMixed` / `RouteMeander` `total_length` is not exact and fillets shorten the path: iterate on the **filleted** centre-line length (`path_length`) until it equals the target.
* Use one fillet radius everywhere and assert every inner straight segment is ≥ 2 × fillet.
* **Resonator meander trace-to-trace spacing (DRC B7):** keep adjacent parallel resonator meander sections (including the coupling leg to the first meander turn) at least **300 µm trench-to-trench** (ground strip width between parallel CPW cuts ≥ 300 µm, i.e. centre-to-centre pitch ≥ `300 µm + W + 2S` = 322 µm for a 10/6 µm CPW).
* **Side-by-side CPW runs:** avoid long parallel runs between different resonators; estimate any remaining parallel run with `purcell_readout_model.parallel_j` (`code/resonator/cpw_pair_cap2d.py`) against target J. `400 µm` remains the default `line_spacing`, `jog`, and `meander_pitch` (not a hard spacing limit). Every CPW still requires ground bonds ([em_simulation.md §3](em_simulation.md#3-ground-wire-bonds-in-every-hfss--q3d-run)).

### 2.3 Verification after a change
* Layout script runs without short-segment warnings; all pin `middle` and `normal` values verified.
* No crossing routes in `fig/chip.png`; target vs. filleted resonator lengths printed equal.
* Check minimum metal width (`code/layout/check_metal_width.py`) and DRC ([chip_drc.md](chip_drc.md)).
* Re-run affected Q3D / HFSS checks ([em_simulation.md](em_simulation.md), [cap_matrix_analysis.md](cap_matrix_analysis.md)) before quoting numbers.

### 2.4 Output file naming
* Write only the GDS files allowed by [gds_export.md §2.3](gds_export.md#23-which-gds-files-may-be-written) (`_sim`, `_fab`, delivery files); no preview, marker-only, or copied GDS.
* GDS and parameter files use `<Design>_<date>` (plus revision/role suffix per [gds_export.md](gds_export.md)), e.g. `gds/<Design>_<date>_wafer_rN_fab.gds`, `gds/<Design>_<date>_wafer_rN_sim.gds`, `parameters/<Design>_<date>.json`.
* Layout figures follow [layout_figures.md](layout_figures.md) (`fig/chip.png`, `fig/<part>.png`).
* Select adopted runs/files explicitly in `design.json`—never infer the adopted result from the newest filename ([parameter_provenance.md](parameter_provenance.md)).
* Remove superseded GDS/figures/logs from active folders (git history preserves them).
