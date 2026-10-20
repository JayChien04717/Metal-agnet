# Qiskit Metal Coding Skills & Best Practices

This document defines the standard coding conventions, parameter parsing guidelines, geometry registration workflows, and pin configuration methods for writing and maintaining custom Qiskit Metal planar components (`QComponent`). Developers and AI agents must strictly follow these rules when implementing or refactoring files such as `pattern design/my_fluxonium.py` or `pattern design/my_tmon.py`.

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
exported GDS, so fillets and junction gaps look exactly as they will on the mask.

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
* Simulation versions keep exactly the same pocket / pad / line geometry and only replace JJs + arrays by one lumped inductor (`fluxonium_forprocess.py`). Flux lines that are shorted to ground are drawn as part of the ground plane there (`flux_as_ground`), to avoid Q3D net-overlap errors.

### 8.3 Routing rules
* `RouteAnchors` with explicit anchor points for control lines (Pathfinder fails easily in a crowded chip); `RouteStraight` between neighbouring tees.
* Readout resonators: `RouteMixed` with one meandered segment below the tee and straight segments for the jog to the qubit; its `total_length` is not exact and fillets shorten the line, so iterate `total_length` on the **filleted** centre-line length (`filleted_length`, `path_length`) until it matches the target.
* 99 µm fillets everywhere; after building, assert every inner path segment is ≥ 2 × fillet.
* Ports: P7/P8 feedline, P1/P6 FX flux, P2/P5 FX charge, P3 TC flux, P4 TC charge; check visually that no routes cross.

### 8.4 GDS export
* `export_gds(design, file)`: no junction cells by default (`include_junctions=False`) so every junction site is an empty gap; arrays only with `include_array=True`.
* `fab=True` applies the lab cheesing settings (`note/gds_export_settings.md`).
* After export, unreferenced cells and dangling references are pruned (`_prune_gds`).

### 8.5 Verification after a change
* `python ftf_layout.py` runs without warnings about short segments.
* Pins / normals printed and sane; no crossing routes in `fig/ftf_layout_chip.png`.
* Resonator lengths: target vs filleted length printed equal.
* Re-run the relevant Q3D / HFSS check (see `skill_lom.md`) before quoting numbers.
