# Qiskit Metal Coding Skills & Best Practices

This document defines the standard coding conventions, parameter parsing guidelines, geometry registration workflows, and pin configuration methods for writing and maintaining custom Qiskit Metal planar components (`QComponent`). Developers and AI agents must strictly follow these rules when implementing or refactoring files such as `pattern design/my_fluxonium.py` or `pattern design/my_tmon.py`.

### Where things live

`skill/` and `note/` hold **general** rules and methods only. Everything about one Design (numbers,
tuning logs, run records, design-specific workflows) lives in that Design's folder.

| folder / file | content |
|---|---|
| `skill/skill.md` | This file: component coding rules, layout-script workflow, environment, simulation rules |
| `skill/gds_export.md` | GDS export: junction-free rule (`jj_metal=False`), cheesing, units, test junctions, checks |
| `skill/layout_figures.md` | Layout figures from the exported GDS: rules, format, adding a view |
| `skill/cap_matrix_analysis.md` | Procedure: Q3D matrix → Ec, J, g, β, κ |
| `skill/mutual_inductance_analysis.md` | Procedure: flux line → loop M from HFSS |
| `skill/chip_drc.md` | Layout DRC before a GDS export |
| `note/` | General methods and shared code: `cap_matrix.py`, `purcell_filter_method.md`, `wallraff reference/` |
| `pattern design/` | Shared QComponents and `gds_tools.py` (all Designs) |
| `Design<N>/notes/` | Design-specific methods, results and tuning logs |
| `Design<N>/parameters/`, `design_parameters.txt` | Final numbers of that Design; never in `note/` or `skill/` |

---

## 🧭 Table of Contents
1. [Custom QComponent Template](#1-custom-qcomponent-template)
2. [Parameter Parsing & Units Evaluation](#2-parameter-parsing--units-evaluation)
3. [Shapely Drawing & Coordinate Transformations](#3-shapely-drawing--coordinate-transformations)
4. [QGeometry Registration (Poly, Path, Junction)](#4-qgeometry-registration-poly-path-junction)
5. [Pin Declaration & Automated Routing](#5-pin-declaration--automated-routing)
6. [Ground Plane Cutout (Subtract Pocket) Implementation](#6-ground-plane-cutout-subtract-pocket-implementation)
7. [Common Pitfalls & Review Checklist](#7-common-pitfalls--review-checklist)
8. [Layout Script Workflow](#8-layout-script-workflow)

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

## 8. Layout Script Workflow

Each Design builds its chip with one script in `Design<N>/pattern design/` (no notebook): a
`build_design()` (or `build_chip()`) that returns the design, `export_gds()` for `gds/`, and
`plot_gds()` for the figures, both from the shared `pattern design/gds_tools.py`. Figures are drawn
from the exported GDS, so fillets and junction gaps look exactly as on the mask.
Design-specific knobs, ports and conventions are in `Design<N>/notes/`.

### 8.1 Knobs
* Put every tunable number in a module-level dict or constant at the top of the script, and let
  `build_design(...)` take overrides, so tuning sweeps never edit the code.
* Derive dependent geometry (e.g. coupler size from qubit pad size) instead of setting it twice.

### 8.2 Component conventions
* Draw in local coordinates, then mirror / rotate / translate with one helper, and create pins through
  a helper that swaps the point order when mirrored, so normals stay outward.
* Check every pin's `middle` and `normal` after a change (`design.components[name].pins`).
* Boolean results (`draw.subtract`) can leave zero-area slivers; drop them before `add_qgeometry`,
  otherwise Ansys fails with *Cannot create rectangle with zero length/width*.
* Junction-table lines: draw all junctions of a device in the same direction (Dolan bridges need parallel junctions).
* **One simulation model per qubit type for every run** (Q3D, HFSS eigenmode, mutual inductance), with the
  same pads / stems / pocket / lines as the mask. Fluxonium: `FluxoniumArrayJJ` (`fluxonium.py`): the JJ is a
  short gap of the lead line at the qubit centre; the array is a U of two junction-table sides joined by a
  metal closing end (an HFSS lumped sheet must end on metal; in Q3D the closing end is a floating net,
  Kron-reduced by `note/cap_matrix.py`). Do not use a single straight inductor.
* Flux lines shorted to ground are drawn as ground for Q3D (`flux_as_ground`), to avoid net-overlap errors.

### 8.3 Routing rules
* `RouteAnchors` with explicit anchor points for control lines (Pathfinder fails easily in a crowded chip);
  `RouteStraight` between neighbouring tees; deterministic serpentines for long resonators.
* `RouteMixed` / `RouteMeander` `total_length` is not exact and fillets shorten the line: iterate on the
  **filleted** centre-line length (`path_length`) until it matches the target.
* One fillet radius everywhere; assert every inner path segment is >= 2 x fillet.
* Every pair of side-by-side CPW runs >= 400 um apart (section 11). Check visually that no routes cross.

### 8.4 GDS export
* Junction-free GDS, cheesing with `fab=True`, unreferenced cells pruned: full rules in `gds_export.md`.

### 8.5 Verification after a change
* The layout script runs without warnings about short segments.
* Pins / normals sane; no crossing routes in the chip figure.
* Resonator lengths: target vs filleted length printed equal.
* Re-run the relevant Q3D / HFSS check (`cap_matrix_analysis.md`) before quoting numbers.

---

## 9. Environment (do NOT search for Python every time)

* **Python for every script in this repo (layout, GDS, Q3D/HFSS drivers):** conda env `qmetal`
  (qiskit_metal, gdspy, geopandas, shapely, pyEPR). Lab PC: `C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe`;
  laptop (user zx123): `C:\Users\zx123\miniforge3\envs\qmetal\python.exe`. AEDT 2021 R2 on both.
  Other envs (`C:\Users\QEL\.venv`, `quantum-metal\.venv`, system python) lack geopandas / gdspy: do not use them.
* Run the layout script from the Design's `pattern design/` folder, e.g. in PowerShell:
  `& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I <layout script>.py`.
* qiskit-metal prints a harmless `Fake_Junctions.GDS` warning on stderr; PowerShell then shows a
  NativeCommandError header, ignore it and read the printed lines below it.

---

## 10. HFSS lessons from past runs (narrow resonances, Purcell channels)

* Driven adaptive on narrow lossless resonances never settles in Delta S (it jumps 0.08-1.8 pass to pass on a
  -60 dB null): do not wait for Delta S; use the fixed passes of section 12.
* Mesh from eigenmode where possible: eigenmode design with the same geometry and ports as 50 Ohm RLC sheets,
  then a driven design in the same project with Mesh Link -> `<eig design> : Setup : LastAdaptive`,
  1 pass. AEDT can refuse a Mesh Link silently (ImportMesh stays false): check it.
* First-order basis for eigenmode; second order ran out of 32 GB at pass 1.
* **Locate first:** an Interpolating sweep can miss a resonance. Run a wide sweep (e.g. 5.8-7.6 GHz, 1 MHz)
  first, then a narrow one (<= 0.25 MHz) around the dips found. Discrete 5 MHz steps miss a filter dip.
* A far-detuned resonator (kHz linewidth) is invisible in a driven sweep: find it with eigenmode.
* Keep >= 80 GB free on C: (HFSS scratch `%TEMP%\hfss_*.pjt` reaches 20-45 GB per session); a full disk
  gives COM error -2147352567 at Analyze. The same COM error also comes from a disconnected net or from
  running out of RAM: read the AEDT message window before guessing.
* All runs of a Design share one AEDT project (`save_aedt.use_project`).

---

## 11. CPW spacing rule: every CPW pair >= 400 um apart (Jay, 2026-11-12)

* Any two CPW runs that lie side by side (parallel, overlapping, not on the same straight line) must be
  >= 400 um centre to centre. This covers the meander pitch, the meander-to-main-line clearance, the
  filter-resonator spacing and the distance to the feed line. Exempt: capacitor electrodes, collinear
  pieces of one line, tees / perpendicular segments, and the hammer stubs (they ARE the coupler C_j).
* Why: at 230-300 um, runs couple to each other and the etch cuts the ground. HFSS then shows couplings
  the capacitor model cannot hold (~1-2 MHz filter-resonator coupling without any coupler).
* Check with a spacing script on the built design (all path pairs, centre to centre) after every layout change.
* Routing that meets the rule: vertical serpentines with runs in 400 um steps, lowest runs >= 650 um above the
  feed line.
* Every CPW still needs ground bonds (section 13; physics in `note/purcell_filter_method.md` section 5).

---

## 12. Fixed pass counts and mesh (Jay, 2026-11-12) — USE THESE, do not wait for convergence

| solver | max passes | initial mesh (length based) |
|---|---|---|
| **Q3D** (capacitance matrix) | **8** | Q3D auto |
| **HFSS eigenmode** | **6** | **10 um** on every CPW line and capacitor electrode |
| **HFSS driven S21** | **3** | **10 um** on every CPW line and capacitor electrode |

* Why: with a 10 um length mesh seeded on the metal from the start, HFSS does not need to refine the mesh
  pass by pass; a few passes are enough, and runs stay short (older 8-20 pass / 5 um runs took hours and filled memory).
* Make these the defaults of every simulation script; deviate only deliberately, and record it in the Design's notes.
* Exception, a full 10 x 10 mm chip with several channels: a 10 um seed mesh does not fit in 24 GB RAM (also not with
  30 um on the lines). Tune per channel or channel pair; run the full chip only as a coarse check.

---

## 13. Ground wire bonds in every HFSS / Q3D run (Jay, 2026-11-12)

* Use qiskit-metal's built-in `hfss_wire_bonds`: set `hfss_wire_bonds = True` on every etched path in the path
  qgeometry table and `sim.renderer._options.update(wb_threshold='350um', wb_offset='0um', wb_size=2)` before `run`.
  The default 400 um skips most of the layout; 350 um gives a bond in each 400 um serpentine bottom and one every
  ~300-350 um on long runs.
* Pieces shorter than the threshold (coupler stubs, tee arms) get extra PEC straps drawn after `render_design`.
* The same bonds in Q3D, eigenmode and driven; Q3D counts the bond nets as ground (a bond net that does not
  touch ground on both sides makes Q3D fail: then run Q3D without bonds, it barely changes series C).
* Plot the bond positions before a long run.

---

## 14. File naming of outputs (Jay, 2026-11-13)

* GDS, figures and parameter files are named `<Design folder>_<date>`, e.g. `gds/DesignN_2026-11-13.gds` (fab, cheesed),
  `gds/DesignN_2026-11-13_nocheese.gds`, `fig/DesignN_2026-11-13_layout.png`, `parameters/DesignN_2026-11-13.json`.
* One parameter file per Design and date in `Design<N>/parameters/`, sections: substrate, chip, cpw, capacitors,
  routing, lengths / channels, wire bonds, solver rules, results, fab. Scripts read the newest one.
* Superseded GDS / figures / logs are removed (git history keeps them).
