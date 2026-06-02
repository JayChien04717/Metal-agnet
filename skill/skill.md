# Qiskit Metal Coding Skills & Best Practices

This document defines the standard coding conventions, parameter parsing guidelines, geometry registration workflows, and pin configuration methods for writing and maintaining custom Qiskit Metal planar components (`QComponent`). Developers and AI agents must strictly follow these rules when implementing or refactoring files such as `fluxonium.py` or `grounded_transmon_coupler.py`.

---

## 🧭 Table of Contents
1. [Custom QComponent Template](#1-custom-qcomponent-template)
2. [Parameter Parsing & Units Evaluation](#2-parameter-parsing--units-evaluation)
3. [Shapely Drawing & Coordinate Transformations](#3-shapely-drawing--coordinate-transformations)
4. [QGeometry Registration (Poly, Path, Junction)](#4-qgeometry-registration-poly-path-junction)
5. [Pin Declaration & Automated Routing](#5-pin-declaration--automated-routing)
6. [Ground Plane Cutout (Subtract Pocket) Implementation](#6-ground-plane-cutout-subtract-pocket-implementation)
7. [Common Pitfalls & Review Checklist](#7-common-pitfalls--review-checklist)
8. [FTF Notebook Generation Workflow](#8-ftf-notebook-generation-workflow)

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

## 8. FTF Notebook Generation Workflow

Use this workflow when generating or modifying Fluxonium-Transmon-Fluxonium (`FTF`) Qiskit Metal notebooks.

### 8.1 Required Clarifying Questions
Before generating a new FTF notebook, ask only for missing or ambiguous values:

1. Chip size and port keep-in distance.
2. Launchpad port assignment (`P1` through `P8`).
3. Left/right fluxonium positions, orientations, flux-line side, and charge-line side.
4. Transmon coupler type and `cross_length_x` / `cross_length_y`.
5. Transmon flux-line port and charge-line port.
6. Readout tee height, waveguide jog points, and resonator length.
7. Whether to generate a preview image and run geometry sanity checks.

### 8.2 Chip And Port Placement Rule
For a 10 mm x 10 mm chip, the chip boundary is `x/y = +/-5.0mm`.

* Launchpads should be placed 0.5 mm inside the chip edge:
  * left edge: `x = -4.5mm`
  * right edge: `x = +4.5mm`
  * bottom edge: `y = -4.5mm`
  * top edge: `y = +4.5mm`
* Edge-following control-route jog points should stay further inside, typically `+/-4.35mm`, to avoid crossing the chip edge after CPW width, gap, fillet, or launchpad geometry is included.
* Update `design.chips.main.size.size_x` and `size_y` to match the real chip size; do not leave old package sizes such as `10.4mm` in the notebook.

### 8.3 Fluxonium Module Rules
* Use `FluxoniumModule` for each fluxonium qubit.
* Keep the FTF-facing pad exposed through the `transmon_coupler` pin.
* Put the flux-bias line on the outer side of each fluxonium by default with `flux_line_options=dict(side="outer", ...)`.
* Place the fluxonium charge line outside the pocket. Compute its pad center from the pocket edge, not from the capacitor pad edge:
  ```python
  pocket_center_y = p.pad_y - p.pad_height / 2 - p.array_leg_length / 2
  pocket_edge_y = pocket_center_y + loc_h * p.pocket_height / 2
  pad_y = pocket_edge_y + loc_h * (opts.pad_sep + opts.pad_height / 2)
  ```

### 8.4 Transmon Coupler Rules
* Prefer a local `TransmonCrossFL`-based coupler subclass for FTF layouts instead of a simplified grounded bar.
* The transmon cross must support independently tunable arm lengths:
  * `cross_length_x`
  * `cross_length_y`
* Do not use a top claw/connector port for the transmon charge or drive line.
* Keep only direct side stubs for coupling to fluxonium pads:
  * `left_couple`
  * `right_couple`
* Add a separate transmon charge/drive line next to the transmon flux line. It should expose `drive` and/or `charge` pins, while the flux line exposes a `flux` pin.

### 8.5 Standard Port Assignment Example
Use this assignment unless the user overrides it:

```python
port_assignment = {
    "P1": "left fluxonium flux bias",
    "P2": "left fluxonium charge",
    "P3": "transmon coupler flux bias",
    "P4": "transmon coupler charge",
    "P5": "right fluxonium charge",
    "P6": "right fluxonium flux bias",
    "P7": "right readout waveguide",
    "P8": "left readout waveguide",
}
```

### 8.6 Readout Waveguide And Tee Rules
* Resonator tee `x` positions do not need to move just to widen the waveguide. Keep the tee x positions near the resonator coupling locations unless the user asks otherwise.
* Tee height can be moved independently, for example:
  ```python
  LR_tee: pos_x="-0.75mm", pos_y="2.5mm"
  RR_tee: pos_x="+0.75mm", pos_y="2.5mm"
  ```
* Use Manhattan jogs in the feed routes to widen the waveguide. Example:
  ```python
  # Left feed: jog from top port to x=-3mm, drop to y=2.5mm, then return to the tee.
  points=[["-1.1mm", "3.0mm"], ["-3.0mm", "3.0mm"], ["-3.0mm", "2.5mm"], ["-0.75mm", "2.5mm"]]

  # Right feed: leave the tee, jog to x=+3mm, rise to y=3mm, then return toward the top port.
  points=[["0.75mm", "2.5mm"], ["3.0mm", "2.5mm"], ["3.0mm", "3.0mm"], ["1.1mm", "3.0mm"]]
  ```

### 8.7 RouteMeander Resonator Rules
When the resonator/cavity needs a controlled entry jog, use `RouteMeander` with `lead.start_jogged_extension`.

```python
from collections import OrderedDict

lr_jogs_in = OrderedDict()
lr_jogs_in[0] = ["R", "200um"]
lr_jogs_in[1] = ["R", "200um"]

rr_jogs_in = OrderedDict()
rr_jogs_in[0] = ["L", "200um"]
rr_jogs_in[1] = ["L", "200um"]

LR_resonator = RouteMeander(
    design,
    "LR_resonator",
    options=dict(
        pin_inputs=dict(
            start_pin=dict(component="LR_tee", pin="second_end"),
            end_pin=dict(component="FL", pin="readout"),
        ),
        total_length=readout_resonator_length,
        trace_width="cpw_width",
        fillet=route_fillet,
        meander=dict(spacing="200um", asymmetry="0um"),
        lead=dict(
            start_straight="1000um",
            end_straight="260um",
            start_jogged_extension=lr_jogs_in,
        ),
    ),
)
```

Mirror the jog directions for the right resonator.

### 8.8 Verification Requirements
After editing or generating an FTF notebook:

* Validate notebook JSON:
  ```powershell
  & 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe' -c "import json; json.load(open('FTF.ipynb', encoding='utf-8')); print('FTF.ipynb JSON OK')"
  ```
* Run `py_compile` on local component files.
* Instantiate the relevant design cells in the `qmetal` environment when possible.
* Print and inspect component pins:
  * `FL`: should include `charge`, `flux`, `readout`.
  * `FR`: should include `charge`, `flux`, `readout`.
  * `TC`: should include `charge`, `drive`, `flux`.
* Check that route names match the intended port assignment.
* Search for stale geometry constants such as `10.4mm`, `5.05mm`, `5.025mm`, or `4.92mm` after converting to a 10 mm chip.
