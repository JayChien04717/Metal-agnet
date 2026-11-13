# Qiskit Metal Coding Skills & Best Practices

This document defines the standard coding conventions, parameter parsing guidelines, geometry registration workflows, and pin configuration methods for writing and maintaining custom Qiskit Metal planar components (`QComponent`). Developers and AI agents must strictly follow these rules when implementing or refactoring files such as `pattern design/my_fluxonium.py` or `pattern design/my_tmon.py`.

### Where things live (SQL HFSS/)

| folder / file | content |
|---|---|
| `skill/skill.md` | This file: component coding rules, FTF layout workflow, environment, HFSS settings |
| `skill/gds_export.md` | GDS export: **junction-free** rule (`jj_metal=False`), cheesing, units, checks |
| `skill/layout_figures.md` | How to make the layout figures ("screenshots") in `Design<N>/fig/`: format, standard set, windows, adding a view |
| `skill/cap_matrix_analysis.md` | Procedure: Q3D matrix → Ec, J, g, β, κ |
| `skill/mutual_inductance_analysis.md` | Procedure: flux line → loop M from HFSS |
| `skill/chip_drc.md` | Layout DRC before a GDS export |
| `note/` | Shared methods with results / tuning logs: `capacitance_method.md` (+ `cap_matrix.py`), `mutual_inductance_method.md`, `zz_method.md`, `purcell_filter_method.md`, `wallraff reference/` |
| `pattern design/` | Shared QComponents (all Designs) |
| `Design<N>/design_parameters.txt` | Final numbers of that Design (Design1 – Design3); never in `note/` |

---

## 🧭 Table of Contents
1. [Custom QComponent Template](#1-custom-qcomponent-template)
2. [Parameter Parsing & Units Evaluation](#2-parameter-parsing--units-evaluation)
3. [Shapely Drawing & Coordinate Transformations](#3-shapely-drawing--coordinate-transformations)
4. [QGeometry Registration (Poly, Path, Junction)](#4-qgeometry-registration-poly-path-junction)
5. [Pin Declaration & Automated Routing](#5-pin-declaration--automated-routing)
6. [Ground Plane Cutout (Subtract Pocket) Implementation](#6-ground-plane-cutout-subtract-pocket-implementation)
7. [Common Pitfalls & Review Checklist](#7-common-pitfalls--review-checklist)
8. [FTF Layout Workflow (ftf_layout.py)](#8-ftf-layout-workflow-ftf_layoutpy)

---

## 1. Custom QComponent Template

All custom planar components must inherit from `QComponent`. A standard `QComponent` implementation consists of three core sections:

```python
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import QComponent

class CustomQubit(QComponent):
    """
    Component Description: Summarize the design purpose and main geometric features of this component.
    """
    
    # 1. Default Options (all tunable geometric dimensions, position, orientation, and layer)
    default_options = Dict(
        pos_x='0um',
        pos_y='0um',
        orientation='0',         # Rotation angle in degrees
        pad_width='120um',
        pad_height='30um',
        connection_gap='6um',
        layer='1',               # Default metal layer
    )
    
    # 2. Component Metadata (defines GUI representations and database mappings)
    component_metadata = Dict(
        short_name='custom_qubit',
        _qgeometry_table_poly='True',
        _qgeometry_table_path='True',
        _qgeometry_table_junction='True' # Set to True if the component features JJs
    )

    # 3. Core Drawing Logic
    def make(self):
        p = self.p  # Essential: Access parsed options as float values in standard SI units (meters)
        
        # Geometry generation and registration logic...
```

---

## 2. Parameter Parsing & Units Evaluation

Qiskit Metal supports flexible unit representations (e.g., `'100um'`, `'2mm'`, `'500nm'`) and string-based mathematical expressions.

### 2.1 🚫 The WRONG Way (Never read `self.options` directly for arithmetic)
```python
# ERROR Example
pad_w = self.options.pad_width  # Retrieves the raw string '120um'
new_w = pad_w * 2               # Python duplicates the string: '120um120um', causing geometry crashes!
```

### 2.2 📐 The RIGHT Way (Use `self.p` or `self.parse_value`)
* **Use `self.p`**: Within `make()`, `self.p` automatically parses and evaluates the string values inside `self.options` into **float values in standard SI units (meters)**.
  ```python
  pad_w = self.p.pad_width   # '120um' is automatically parsed as 0.00012 (float)
  new_w = pad_w * 2          # Safe arithmetic operations (results in 0.00024)
  ```
* **Manual Parsing of Dynamic Values**: If you need to evaluate an arbitrary string expression during execution, use `self.parse_value()`:
  ```python
  custom_gap = self.parse_value('6um')  # Returns 6e-6 (float in meters)
  ```

---

## 3. Shapely Drawing & Coordinate Transformations

Planar components in Qiskit Metal are built using **Shapely geometry objects**.

### 3.1 Common Shapely Geometries
* **Rectangle**:
  ```python
  # Parameters: width, height, center_x, center_y
  rect = draw.rectangle(p.pad_width, p.pad_height, 0, 0)
  ```
* **Rounded Rectangle (Fillet)**:
  ```python
  rect = draw.rectangle(p.pad_width, p.pad_height, 0, 0)
  # Apply buffer padding to round corners (join_style=1 stands for round)
  rounded_rect = rect.buffer(-p.connection_gap, join_style=1).buffer(p.connection_gap, join_style=1)
  ```
* **Polygon**:
  ```python
  poly = draw.Polygon([(x1, y1), (x2, y2), (x3, y3)])
  ```
* **LineString (Path)**:
  ```python
  line = draw.LineString([(x1, y1), (x2, y2)])
  ```

### 3.2 🛠️ Standard Coordinate Transformation Flow (Rotate then Translate)
To ensure the component responds correctly to the user's `pos_x`, `pos_y`, and `orientation` parameters, all drawn Shapely geometries **must be rotated and then translated** prior to registration:

```python
# 1. Draw geometry in the local coordinate system (centered at 0, 0)
geom = draw.rectangle(p.pad_width, p.pad_height, 0, 0)

# 2. Perform transformations (ROTATE first, then TRANSLATE)
geom = draw.rotate(geom, p.orientation, origin=(0, 0)) # orientation is in degrees
geom = draw.translate(geom, p.pos_x, p.pos_y)          # Translate to the actual chip position

# 3. Register geometry
self.add_qgeometry('poly', 'my_pad', geom, layer=p.layer)
```

---

## 4. QGeometry Registration (Poly, Path, Junction)

Once drawn and transformed, all geometry objects must be registered via `self.add_qgeometry` to appear in GUI rendering or layout exports.

| Kind | Target | Shapely Type | Mandatory Parameter & Description |
| :--- | :--- | :--- | :--- |
| **`'poly'`** | Capacitor pads, ground cutouts, launch pads | `Polygon` / `MultiPolygon` | Standard surface geometry. Use `subtract=True` for ground plane cutouts. |
| **`'path'`** | CPW lines, wire bonds, control lines | `LineString` / `MultiLineString` | Requires `width` parameter (usually `self.p.trace_width`) to define the trace thickness. |
| **`'junction'`** | Josephson junctions, SQUID loops | `LineString` (typically a very short bridge) | Requires `width` parameter. Evaluated as a lumped element boundary in EM simulations. |

### 4.1 Implementation Examples
```python
# Register solid metal capacitor pad (Poly)
self.add_qgeometry('poly', 'left_pad', left_pad_geom, layer=p.layer)

# Register CPW feedline (Path)
self.add_qgeometry('path', 'feedline', line_geom, width=p.trace_width, layer=p.layer)

# Register Josephson Junction (Junction)
# Drawn as a tiny line segment bridging the gap between capacitor pads
jj_line = draw.LineString([(0, -gap/2), (0, gap/2)])
jj_line = draw.rotate(jj_line, p.orientation, origin=(0, 0))
jj_line = draw.translate(jj_line, p.pos_x, p.pos_y)
self.add_qgeometry('junction', 'rect_jj', jj_line, width=p.jj_width, layer=p.layer)
```

---

## 5. Pin Declaration & Automated Routing

Declaring Pins allows routing algorithms (e.g., `BaseRoute` or meandered line classes) to establish connections between components.

### 5.1 Three Pillars of a Pin
1. **Location & Width (LineString)**: A line segment perpendicular to the port's exit direction (representing the trace width at the boundary).
2. **Outward Normal Vector (Normal)**: A 2D vector pointing **outwards** from the component, defining the initial route direction.
3. **Trace Width & Gap (Width & Gap)**: Propagated CPW dimensions passed to the automatic routing engine.

### 5.2 Standard Pin Declaration Example
The following code defines and registers a Pin on the right edge of a component facing outwards (`normal = (1, 0)`):

```python
# 1. Define Pin center and boundary line in the local coordinate system
# Center located at (pad_w / 2, 0)
pin_center_x = pad_w / 2
pin_center_y = 0

# Define a short LineString representing the trace width
pin_line = draw.LineString([
    (pin_center_x, pin_center_y - p.trace_width/2), 
    (pin_center_x, pin_center_y + p.trace_width/2)
])

# 2. Define outward normal vector (facing positive X-axis)
normal_vector = draw.Point(1, 0)

# 3. Apply standard rotation and translation transformations
pin_line = draw.rotate(pin_line, p.orientation, origin=(0, 0))
pin_line = draw.translate(pin_line, p.pos_x, p.pos_y)

# Normal vector only requires rotation
normal_rotated = draw.rotate(normal_vector, p.orientation, origin=(0, 0))
normal_coords = list(normal_rotated.coords)[0] 

# 4. Register Pin
self.add_pin(
    name='readout_port',
    points=pin_line,
    addressable_points=pin_line.coords,
    normal=normal_coords,              # Pass rotated normal coordinates (x, y)
    width=p.trace_width,
    gap=p.trace_gap,
    input_as_norm=True                 # Confirms that the normal parameter represents a vector
)
```

---

## 6. Ground Plane Cutout (Subtract Pocket) Implementation

In a single-layer planar chip design, ground plane cutouts (pockets) must be generated around capacitor pads to prevent electrical shorts.

* **Pocket Dimensions**: The cutout area should be larger than the physical metal pads by at least `connection_gap`.
* **Subtract Overlays**:
  * Registering a Poly geometry with `subtract=True` instructs Qiskit Metal to remove ground plane copper inside that boundary.
  * Physical metal shapes (registered with `subtract=False` or default) placed inside the cutout boundary are preserved, creating a floating metal structure surrounded by an etched vacuum gap.

```python
# 1. Draw physical capacitor pad
pad = draw.rectangle(p.pad_width, p.pad_height, 0, 0)

# 2. Draw ground cutout pocket (larger than pad)
pocket = draw.rectangle(p.pad_width + 2*p.connection_gap, p.pad_height + 2*p.connection_gap, 0, 0)

# 3. Transformations
pad = draw.translate(draw.rotate(pad, p.orientation), p.pos_x, p.pos_y)
pocket = draw.translate(draw.rotate(pocket, p.orientation), p.pos_x, p.pos_y)

# 4. Registration
self.add_qgeometry('poly', 'qubit_pad', pad, layer=p.layer)  # Physical metal (kept)
self.add_qgeometry('poly', 'qubit_pocket', pocket, subtract=True, layer=p.layer) # Etched ground pocket (cutout)
```

---

## 7. Common Pitfalls & Review Checklist

### 🚨 Common Pitfalls
1. **Josephson Junction Shorts & Open Circuits**:
   * JJs are represented as **LineStrings** in the `junction` table.
   * Make sure both endpoints of the LineString **overlap or touch** the solid metal capacitor pads. Gaps between the junction line and the pads will prevent proper port excitation in EM simulations.
2. **Self-Intersecting Polygons**:
   * Hand-calculated coordinates for custom shapes (e.g. U-shaped paths, asymmetric shunt pads) must never cross or fold back on themselves. Shapely will throw `TopologyException` errors during parsing.
3. **Invalid Outward Normals**:
   * Normal vectors must be non-zero (typically unit vectors). If normal coordinates are pointing inwards, automatic route classes will draw CPW routes directly through the qubit body.

### 📋 Code Review Checklist
- [ ] Every geometric calculation in `make()` utilizes `self.p` (floats in meters) instead of the raw `self.options`.
- [ ] All Shapely shapes are rotated first (`draw.rotate`) and then translated (`draw.translate`) before registration.
- [ ] Cutout pockets are registered with `subtract=True`, whereas solid metal shapes retain the default value.
- [ ] All Pin normals point outwards and are properly rotated according to the component `orientation`.
- [ ] Josephson junctions are registered in the `junction` table rather than the standard `path` table.

---

## 8. FTF Layout Workflow (ftf_layout.py)

The chip is built by `pattern design/ftf_layout.py` (no notebook). `build_design()`
returns the design; `export_gds()` writes `gds/`; `plot_gds()` draws figures from the
exported GDS (both in the shared `pattern design/gds_tools.py`, see `gds_export.md`), so fillets and junction gaps look exactly as they will on the mask.

### 8.1 Where the knobs are
* `FX_OPTIONS` (shared fluxonium geometry) and `FX_PER_QUBIT` (per-qubit overrides, e.g. the readout gap that sets g for each resonator frequency).
* `TMON_OPTIONS`: coupler. `pad_arm_length` and the coupler height are **derived** from the FX pads (`PADDLE_SEP` = paddle-to-pad gap), so changing the FX pad size moves the coupler with it.
* `RESONATORS`: per resonator `f`, `kappa`, `C_pad` (Q3D readout-pad capacitance), `C_extra` (full-circuit Q3D coupling not in the tee fit) and `hfss_scale` (f_HFSS / f_target). The tee coupling length and the total length are computed from these.
* `FX_X`, `Y0`, `FEED_Y`: qubit positions and feedline height.
* `build_design(for_process=True, fx_options=..., tmon_options=..., fx_x=..., fx1_options=..., fx2_options=...)` overrides any of these for tuning sweeps.

### 8.2 Component conventions used here
* Draw in local coordinates, then mirror / rotate / translate with one helper (`_place`), and create pins through a helper that swaps the point order when mirrored, so normals stay outward.
* Check every pin's `middle` and `normal` after a change (`design.components[name].pins`).
* Boolean results (`draw.subtract`) can leave zero-area slivers; drop them before `add_qgeometry` (`_clean` in `my_tmon.py`), otherwise Ansys fails with *Cannot create rectangle with zero length/width*.
* Junction-table lines: draw all junctions of a device in the same direction (Dolan bridges need parallel junctions).
* **One fluxonium simulation model for every run** (Q3D, HFSS eigenmode, mutual inductance): `FluxoniumArrayJJ` (`fluxonium.py`). Same pads / stems / pocket / lines as the mask; the loop is `MyFluxonium._loop`: lead line from the stem ends with the JJ as its short middle gap at the qubit centre (`L_j`, or the port in the mutual run); U = two junction-table sides `array_top` / `array_bot` (`l_inductance / 2` each) joined by a metal closing end `array_end` (an HFSS lumped sheet is straight and must end on metal; the closing end is a floating net in Q3D, Kron-reduced by `cap_matrix.py`). No separate JJ branch. Do not make a fluxonium simulation version with a single straight inductor. Design2 (`FluxoniumPocket`) uses the same model through `u_array='True'` (`build_design(for_process=True)`). Flux lines shorted to ground are drawn as ground for Q3D (`flux_as_ground`), to avoid net-overlap errors.

### 8.3 Routing rules
* `RouteAnchors` with explicit anchor points for control lines (Pathfinder fails easily in a crowded chip); `RouteStraight` between neighbouring tees.
* Readout resonators: `RouteMixed` with one meandered segment below the tee and straight segments for the jog to the qubit; its `total_length` is not exact and fillets shorten the line, so iterate `total_length` on the **filleted** centre-line length (`filleted_length`, `path_length`) until it matches the target.
* 99 µm fillets everywhere; after building, assert every inner path segment is ≥ 2 × fillet.
* Ports: P7/P8 feedline, P1/P6 FX flux, P2/P5 FX charge, P3 TC flux, P4 TC charge; check visually that no routes cross.

### 8.4 GDS export (full rules: `gds_export.md`)
* The GDS is **junction-free**: `export_gds(build_design(jj_metal=False), file)`. `jj_metal='False'`
  drops leads / JJ stubs / array end (FX) and the SQUID block (coupler); `export_gds` empties the
  junction table. Simulations keep the default `jj_metal='True'`.
* A new component with junction-related metal must get a `jj_metal` option too.
* `fab=True` applies the lab cheesing settings; unreferenced cells are pruned (`_prune_gds`).

### 8.5 Verification after a change
* `python ftf_layout.py` runs without warnings about short segments.
* Pins / normals printed and sane; no crossing routes in `fig/layout/ftf_layout_chip.png`.
* Resonator lengths: target vs filleted length printed equal.
* Re-run the relevant Q3D / HFSS check (see `cap_matrix_analysis.md`) before quoting numbers.

---

## 9. Environment (do NOT search for Python every time)

* **Python for every script in this repo (layout, GDS, Q3D/HFSS drivers):** conda env `qmetal`
  (qiskit_metal, gdspy, geopandas, shapely, pyEPR). Lab PC: `C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe`;
  laptop (user zx123): `C:\Users\zx123\miniforge3\envs\qmetal\python.exe`. AEDT 2021 R2 on both.
  Other envs (`C:\Users\QEL\.venv`, `quantum-metal\.venv`, system python) lack geopandas / gdspy: do not use them.
* Run from the `pattern design/` folder, e.g. in PowerShell:
  `& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I ftf_layout.py` (or `purcell_layout.py` for Design4).
* qiskit-metal prints a harmless `Fake_Junctions.GDS` warning on stderr; PowerShell then shows a
  NativeCommandError header, ignore it and read the printed lines below it.
* Design folders: Design1 / Design3 = FTF chip, Design2 = other layout, **Design4 = transmission Purcell filter**
  (`pattern design/purcell_layout.py`, component `interdigital_cap.py`, 6.7 GHz).

---

## 10. HFSS run settings from past runs (Design6, single Purcell channel, CPW 10/6, 2021 R2)

**Rule for every Purcell-filter project (Jay, 2026-11-09): mesh from EIGENMODE, driven does NOT refine.**
1. Eigenmode design with the same geometry (ports drawn as 50 Ω RLC sheets = the same objects as the driven lumped
   ports, same wire bonds): first-order basis (second order ran out of 32 GB at pass 1), initial length-based mesh
   10 µm on the CPW lines / 5 µm on the capacitor electrodes, Delta F 0.1 %, n_modes = number of resonances.
2. Driven design in the SAME AEDT project: Mesh Link → `<eig design> : Setup : LastAdaptive`, Maximum Passes = 1
   (solve once on that mesh), second order, then the sweep (Fast, 0.25 MHz).
3. Script: `Design6/simulation/driven/hfss6_meshlink.py` (`SKIP_EIG=1`, `EIG_DESIGN=...` reuse an eigenmode design
   already solved). All runs of a Design share one AEDT project (`save_aedt.use_project`).

Why: driven adaptive on narrow lossless resonances never settles in ΔS (table below: 15–18 passes, ΔS 0.01–0.08)
and 5 µm initial mesh + driven passes filled the 32 GB.

Older driven-adaptive reference runs:

Reference values measured on the Design6 c6_T3 geometry (filter dip ~6.62 GHz). Use them as the default
instead of re-testing; update this table when a new run changes the conclusion.

| run | type / sweep | step | passes (elements) | last ΔS / ΔF | time | result |
|---|---|---|---|---|---|---|
| c6_T3  | driven, Interp. 6.3–7.1, sol 6.7 GHz | 0.5 MHz | 15 | – | 39 min | dip 6.6230 |
| c6_T3d | driven, Discrete 6.58–6.80 | 5 MHz | 15 | – | 156 min | 6.590 (step too coarse, dip missed) |
| c6_T3g | driven, Discrete 6.56–6.62, sol 6.59 | 1 MHz | 15 (274 k) | 0.011 | – | dip 6.6190 |
| c6_T3h | driven, Interp. 6.45–6.80, sol 6.62 | **0.25 MHz (smallest)** | 18 (604 k) | 0.079 | 137 min | dip 6.6182 |
| c6_T3i | as T3h + 12 um cap mesh | 0.25 MHz | crashed (memory) after 2 h 40 | – | – | – |
| e6_T3b | eigenmode, MINF 6.0 | – | 15 (171 k) | 0.04 % | ~18 min | 6.3419 / 6.6670 |

Conclusions
* **Minimum passes, driven: 15** (~275 k elements). 15 → 18 passes (275 k → 604 k) moved the dip only 0.8 MHz
  (6.6190 → 6.6182) but cost ~100 extra min. ΔS never drops below 0.003 when the solution frequency sits on a
  lossless -60 dB null (ΔS jumps 0.08–1.8 pass to pass), so do not wait for ΔS: fix MAXP=15, min_passes 6.
* **Minimum passes, eigenmode: 14** (~130 k elements): pass 14 is within 1 MHz of pass 15 for both modes
  (pass 13 is still 5–7 MHz off). max_passes 15 is enough.
* **Frequency step:** Interpolating 0.25 MHz over ≤ 350 MHz is the finest used and resolves the filter dip;
  Discrete 5 MHz misses it. For the two-dip (ω_r = ω_p, split ~2J ≈ 17–24 MHz) shape use ≤ 0.25 MHz.
  A far-detuned resonator has a kHz linewidth and is not visible in any of these sweeps: find it with eigenmode.
* Do not add CAP_MESH=12um on top of 18 passes (out of memory). Keep ≥ 80 GB free on C: (HFSS scratch
  `%TEMP%\hfss_*.pjt` reaches 20–45 GB per session; a full C: gives COM error -2147352567 at Analyze).
* **Interpolating sweep can miss a resonance (c6_T5, 2026-11-08):** 6.60–6.80 GHz, 0.25 MHz, 15 passes, solve 6.70 came
  out smooth (−17 dB), while a 5.8–7.6 GHz sweep of the same geometry (12 passes) shows a −51 dB, 32 MHz-wide notch at
  6.754 GHz. Always run a wide locate sweep first (5.8–7.6 GHz, 1 MHz, MAXP 12, ~20 min), then a narrow sweep around
  the dips found; never trust a flat narrow Interpolating sweep on its own.

---

## 11. CPW spacing rule: every CPW pair >= 400 um apart (Jay, 2026-11-12)

* Any two CPW runs that lie side by side (parallel, overlapping, not on the same straight line) must be
  >= 400 um centre to centre. This covers the meander pitch, the meander-to-main-line clearance, the
  filter-resonator spacing and the distance to the feed line. Exempt: capacitor electrodes, collinear
  pieces of one line, tees / perpendicular segments, and the hammer stubs (they ARE the coupler C_j).
* Why: at 230-300 um, runs couple to each other and the etch cuts the ground. HFSS then shows couplings
  the capacitor model cannot hold (Design4: J_par ~1-2 MHz without a hammer).
* Check: `Design4/pattern design/cpw_spacing_check.py` -> `report(design, min_um=400)`. Run it after every
  layout change (`python cpw_spacing_check.py` compares the old and the wide Purcell layout).
* Routing that meets the rule: `purcell_6500.serpentine_tail` (vertical runs in 400 um steps, bottoms >= 650 um
  above the feed line). `build(..., wide=True)` is the default; `wide=False` gives the old 230 um folded tails.
  HFSS / Q3D design names start with `W` for the wide layout (`hfss_hammer_straps.prefix`).
* Every CPW still needs ground straps (`ground_straps(design, 'cpw')`); see note/purcell_filter_method.md.

---

## 12. Fixed pass counts and mesh (Jay, 2026-11-12) — USE THESE, do not wait for convergence

| solver | max passes | initial mesh (length based) | script default |
|---|---|---|---|
| **Q3D** (capacitance matrix) | **8** | Q3D auto | `q3d/q3d_hammer_gap.py --maxp 8` |
| **HFSS eigenmode** | **6** | **10 um** on every CPW line and capacitor electrode | `eigenmode/hfss_hammer_straps.py --maxp 6`, `MESH_LINE = MESH_CAP = '10um'` |
| **HFSS driven S21** | **3** | **10 um** on every CPW line and capacitor electrode | `driven/hfss_hammer_driven.py --maxp 3 --mesh-line 10um --mesh-cap 10um` |

* Why: with a 10 um length mesh seeded on the metal from the start, HFSS does not need to refine the mesh
  pass by pass; a few passes are enough, and runs stay short (older 8-20 pass / 5 um runs took hours and filled memory).
* All three are the script defaults now; pass `--maxp` only to deliberately deviate, and record it in the note.
* Exception, full 10 x 10 mm chip (Design5, 6 channels): 10 um on every line does not fit in 24 GB RAM; even 10 um caps +
  30 um lines ran out of memory at pass 2 (1.2 M elements at pass 1, 2026-11-12). Tune per channel / channel pair instead.
* These supersede the older pass numbers in section 10 (15 driven / 14 eigenmode), which are kept only as history.

---

## 13. Ground wire bonds in every HFSS / Q3D run (Jay, 2026-11-12)

* Use qiskit-metal's built-in `hfss_wire_bonds`: `purcell_6500.set_wire_bonds(design)` flags every etched path;
  `sim.renderer._options.update(P.WB_OPTIONS)` (wb_threshold 350 um, wb_size 3) before `run`. Default 400 um
  skips most of the layout; 350 um gives a bond in the middle of each 400 um serpentine bottom and ~every 300-350 um on long runs.
* Pieces shorter than the threshold (hammer stubs, tee arms) get `draw_hammer_straps(renderer, d, 'all')`
  inside `_analyze` after `render_design`.
* Check the positions before a run: `python plot_wirebonds.py --gap <g>` -> `fig/purcell_wide400_gap<g>_wirebonds.png`.

---

## 14. File naming of outputs (Jay, 2026-11-13)

* GDS, figures and parameter files are named `<Design folder>_<date>`: `gds/Design5_2026-11-13.gds` (fab, cheesed),
  `gds/Design5_2026-11-13_nocheese.gds`, `fig/Design5_2026-11-13_layout.png`, `parameters/Design5_2026-11-13.json`.
* One parameter file per Design and date in `Design<N>/parameters/`, sections: substrate, chip, cpw, capacitors,
  routing, lengths / channels, wire bonds, solver rules, results, fab. Scripts read the newest one.
* Superseded GDS / figures / logs are removed (git history keeps them).
