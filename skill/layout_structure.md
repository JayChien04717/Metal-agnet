# Chip Layout Structure: Placement, Routing, and Build Entry Point

Defines the standard structure of `Design<N>/layout/` so layouts can be modified, re-routed, and adapted across e-beam sites without duplicating code or mutating finished designs.

Related rules: [coding & stage gates](skill.md), [shared QComponents](qcomponents.md), [layout DRC](chip_drc.md), [single-die GDS export](gds_export.md), [wafer deliveries](fabrication_output.md).

---

## 1. Three-file `Design<N>/layout/` structure

Keep component placement, routing, and build logic separate:

```text
Design<N>/layout/
  placement.json   # Fixed components: class, position/orientation, options (incl. site markers, coupons, die ID)
  routes.json      # CPW routes: start/end pins, anchors, fillet, target filleted length, meander rules
  build.py         # Single entry point: build_design(variant='sim'|'fab', overrides=None)
  components/      # Optional Design-specific QComponents only
```

* **`placement.json`** declares every placed component—qubits, couplers, tees, launch pads, site alignment markers, coupons, and die ID—by class name, coordinates, and options.
  * Place the target e-beam site's alignment marker component directly here (e.g. `TSRIA12Mark` or `TiAuChipMarker` from `pattern design/marker/`), with no extra wrapper layer. Changing sites updates or removes this entry in a new layout revision.
  * Because coordinates live in JSON, every layout revision is directly diffable.
* **`routes.json`** declares every CPW route by start/end `(component, pin)`, explicit anchors, fillet radius, target filleted length (for resonators), and meander/spacing rules. Re-routing around a marker or pad edits anchors here while the length solver restores target resonator lengths.
* **`build.py`** reads `placement.json` and `routes.json`, instantiates the components and routes, and returns the Qiskit Metal design.

---

## 2. `build.py` interface and dependency rules

1. **Single entry point:** expose `build_design(variant='sim' | 'fab', overrides=None)` (accepting optional parameter overrides so sweeps never edit source files).
2. **One-way dependency:** simulation scripts (`simulation/**`) and GDS exporters call `layout/build.py`. Layout scripts must **never** import or call simulation modules.
3. **No post-hoc mutation:** instantiate the target classes and routes directly from the JSON definitions. Do **not** monkey-patch or swap component classes, delete components to redraw them, or patch routes after building.
4. **Clean `sim` vs `fab` variants:**
   * `variant='sim'`: solver geometry (`jj_metal=True`, no Dolan contact hooks, no e-beam junction layers, no cheesing). Site markers may still be present in Stage 1 so `fig/chip.png` shows marker openings and passes keep-out checks; exclude non-electrical marker artwork from EM solvers as needed.
   * `variant='fab'`: mask geometry (`jj_metal=False`, `fab_contacts=True`, plus requested coupons, markers, and e-beam junction layers per [gds_export.md §2.1](gds_export.md#21-fab-gds-workflow)).
5. **Revisions:** any layout or site-marker change creates a new revision recorded in `design.json`. Do not store layout scripts inside `Fabrication output/<YYYYMMDD>/`.

---

## 3. Shared layout utilities (`code/general/`)

Two shared utilities belong in `code/general/` rather than being duplicated per Design:

1. **QGeometry & GDS keep-out check (`code/layout/check_keepouts.py`):**
   * `check_keepouts(design, min_um=100)` collects keep-out polygons from marker and coupon components ([qcomponents.md](qcomponents.md)) and measures clearance to all active circuitry and CPW route etches directly in QGeometry during **Stage 1** (before GDS export).
   * CLI for any exported single-die GDS: `python -I -X utf8 -B code/layout/check_keepouts.py <mask.gds> [--min-um 100] [--sites NTHU,TSRI]` (uses `check_gds_keepouts()`, exits non-zero on any violation).
   * Default marker clearance (universal across all Designs): **≥ 100 µm** from **both** NTHU TiAu (`TiAuChipMarker`) and TSRI A12 (`TSRIA12Mark`) bare-rectangle openings to any active route, pad, coupon, or circuit etch (`include_site_markers=('NTHU', 'TSRI')` by default, [fabrication_output.md §2](fabrication_output.md#2-alignment-markers-by-e-beam-site)), so every chip layout is site-portable without re-routing.
2. **Route length compensation solver (`code/layout/route_length_solver.py`, `solve_route_length(...)`):**
   * Given a route's anchor path and target filleted centre-line length, iterates meander / jog geometry until the filleted `path_length` equals the target while enforcing minimum segment length (`≥ 2 × R`) and meander trench-to-trench spacing (`≥ 300 µm`, [chip_drc.md B2, B7, B8](chip_drc.md)).

---

## 4. Shared fab export pipeline

Single-die fab export follows a fixed, design-agnostic sequence:

```text
build_design(variant='fab') → export_gds → cheese (wafer_cheese_grid.py) → verify
```

* Shared export scripts under `code/general/` must **never** hardcode Design names, Design-specific branches, or paths to removed folders.
* All Design-specific options (launch assignments, coupon styles, marker selection, fab hooks) live inside each Design's `layout/` (`placement.json`, `routes.json`, `build.py`) and `design.json`.
