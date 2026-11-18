# Skill: GDS export

How a mask GDS is written for any Design. Design-specific export options and checks are in
`Design<N>/notes/`.

---

## 1. Electrical geometry and fabrication

The fab draws every Josephson junction, junction array and test junction. The mask therefore holds
**no junction cell and no junction-related metal** (leads, JJ stubs, array closing ends, SQUID
blocks), unless the Design explicitly keeps landing electrodes for the fab.

* Every component with junction-related metal has a `jj_metal` option (default `'True'` for
  simulations). The GDS is built from a design made with `jj_metal=False`.
* `export_gds` empties the `junction` qgeometry table, so no JJ library is needed.
* A new component with junction-related metal must get a `jj_metal` option too.
* Q3D and HFSS render the unnotched simulation QComponent directly from QGeometry. The junction
  table provides the lumped elements. A layout export cannot encode those boundary conditions.
* A dedicated simulation GDS, as used in Design3 and Design6 `*_sim.gds`, preserves the solver geometry on layer 1
  and annotates the lumped junction sheets on layer 50/0. In Designs 1, 2, 3 and 3_2 the two GDS
  deliverables are the unnotched base-metal layout and the fabrication mask. The simulation builder
  provides their lumped geometry.
* Use the fabrication QComponent to add contact hooks to both sides of each device junction.
  The current hook is 2 µm wide and 4 µm deep. The two TC ground-side cuts use the same signed
  +2 µm offset in local x. Confirm clear junction gaps and that the hook does not sever the landing.
* Resonator-only Designs 4 and 5 have no device junctions and need no junction hooks.

## 2. How to export

1. Before exporting, ask whether the four corner test-junction structures should be included,
   unless the user already chose for this export. Omitted: remove both the pads and the etched
   pockets. Included: export pads and pockets; the junction itself is left to the fab.
2. Use only the shared code in `code/general/gds_tools.py`: `export_gds`, `apply_fab_gds_options`,
   `plot_gds`, `junction_window`. Never copy these into a Design.
3. `export_gds` edits the design in place (junction table emptied): build a new design for each file.
4. Write a pair of files: the unnotched reference or dedicated simulation GDS, and the fab mask
   with hooks and cheesing (`fab=True`). Use explicit `_sim.gds` and `_fab.gds` names for new
   dedicated simulation/fab exports.

## 3. Checks after an export

1. Open the GDS in KLayout: units µm, chip outline correct (10 × 10 mm chip: −5000 … +5000 µm).
2. Junction sites are empty gaps (or hold only the landing electrodes the Design asks for).
3. Re-export after any component change; commit GDS and figures together.
4. Probe that every etched feature survived the export (e.g. `gdstk` `polygon.contain` at a point inside each
   flux-line / CPW gap). Known qiskit-metal pitfall (Design6, 2026-11-14): a die-ID label that starts with a
   component name (die ID `TCZ00`, coupler `TC`) made the GDS export silently drop that component's flux-line
   and route etches. Never start a die ID with a component name.

Wafer placement exports put the die boundary and device metal on layer 1/0.
The grounded metal perimeter is a 10 mm square, with outer corners at
(±5000, ±5000) µm. Simulation junction annotations retain their separate layer.

For wafer fab exports, the shared [`wafer_cheese_grid.py`](../code/general/wafer_cheese_grid.py)
rebuilds the perforations in the merged positive-metal ground and checks the GDS
in KLayout. The whole hole must fit inside x, y = ±4,452.5 µm, the measured
outermost hole edge of the Design3 reference. A 10 mm die therefore has a
547.5 µm unperforated margin on all four sides. The unnotched and simulation
exports omit this fab operation.

Near a wafer launch, use the local clearance profile measured from the eight
launches of the [Design3 r5 mask](../Design3/gds/previous/Design3_3_2026-11-16_wafer_r5_fab.gds).
Measure the longitudinal distance from each nominal PCB anchor inward to the
hole centre and the transverse distance from the launch centreline to the closest hole edge. The
minimum clearances are the same on both sides of the port.

| Inward distance from PCB anchor (µm) | Minimum transverse clearance to a hole edge (µm) |
|---:|---:|
| 50 | 397.5 |
| 100 | 297.5 |
| 150 | 247.5 |
| 200 | 197.5 |
| 250 | 97.5 |
| 300 | 97.5 |

Interpolate between the measured rows and carry the 50 µm row outward to the
anchor. Apply the profile around each populated launch and apply the four-edge
margin to the entire die. Put 5 × 5 µm holes on a 20 µm centre-to-centre grid,
with a 10 µm coordinate offset on each axis to align its outer rows to the
reference. Keep at least 50 µm from the hole edge to other etches. The grid
gives a 15 µm minimum metal strip between adjacent holes. Use the mask's
existing square or rounded hole style. Standalone test-junction masks use the
same size and pitch and their 500 µm edge keep-out.

Standalone design masks retain their existing design identifier in the upper-left
margin, 50 µm inside both edges. Assembled wafers add the placement-specific IDs
described below.

### 4-inch wafer placement IDs

Every populated die in a 4-inch wafer export carries a unique ID in its upper-right
corner. Order the occupied positions by descending y, then ascending x, reading
each row left to right and the rows top to bottom. Skip empty wafer-edge slots.
Number consecutively from D00 with at least two decimal digits. A 52-die wafer
therefore ends at D51.

Use polygon text on the process layer, with a 100 µm nominal font size and the
glyph's right and top bounds 50 µm inside the actual source die edges. Etch the
text into a solid ground margin. Where the corresponding vendor-die corner is
blank, use positive metal text. Check the whole label window for existing features
before editing and preserve the wafer's street width and alignment crosses.

Record each die's ID, cell, region, centre and local glyph bounds in the assembly
JSON. Keep placement IDs distinct from the source design identifier. The assembler
uses shared source cells and one wrapper per numbered die to keep the repeated
electrical geometry reproducible.

Use `code/general/assemble_wafer.py` for the positive-metal GDS and numbered map.
The optional `--output-stem` selects the companion GDS/JSON/PNG names. Generate
the matching etch GDS/DXF with `code/general/invert_wafer_tone.py --source-stem`
and the required exposure scope. Read both formats back, check the D00-to-last
sequence and upper-right positions, and verify the text polarity after inversion.

Every delivered wafer etch GDS and DXF includes the same 100 µm dicing street
lines. Extend the horizontal and vertical grid lines to the wafer boundary and
exclude the alignment-cross polygons from exposure so that the crosses remain
metal. Adding or updating die IDs must preserve this street geometry. The
inverter defaults to `--dies-and-streets`. Its other scopes reproduce historical
or explicitly requested variants. Compare the street geometry of the final GDS
with the approved street reference and verify the DXF against that GDS.

## 4. Fab / cheesing settings (`fab=True`)

In the code: `apply_fab_gds_options()` in `gds_tools.py`.

```python
gds_render.options.short_segments_to_not_fillet = False
gds_render.options.tolerance = 1e-6
gds_render.options.fabricate = True
gds_render.options.negative_mask = dict(main=[2])
gds_render.options.no_cheese['view_in_file']['main'] = {1: False, 2: False}
gds_render.options.no_cheese['buffer'] = '50um'
gds_render.options.cheese['edge_nocheese'] = '0.5mm'
gds_render.options.cheese['view_in_file']['main'] = {1: True, 2: False}
gds_render.options.cheese.update(dict(cheese_0_x=0.005, cheese_0_y=0.005, delta_x=0.020, delta_y=0.020))
```

Design units are mm, so 0.005 = 5 µm and 0.020 = 20 µm.

| option | value | meaning |
|---|---|---|
| `short_segments_to_not_fillet` | False | Fillet every corner, even on short segments (keep every route segment ≥ 2 × fillet). |
| `tolerance` | 1e-6 m = 1 µm | Tolerance for the filleted path outlines. |
| `fabricate` | True | Write only the final mask; drop intermediate cells (keep-out, hole cell, uncheesed ground). |
| `negative_mask` | `{main: [2]}` | Layer 2 would be negative. Layer 1 (all metal here) stays **positive**. |
| `no_cheese.buffer` | 50 µm | No holes within 50 µm of any etched feature. |
| `no_cheese.view_in_file` | `{1: False, 2: False}` | Keep-out region not written to the file. |
| `cheese.view_in_file` | `{1: True, 2: False}` | Cheese the ground on layer 1 only. |
| `cheese_0_x`, `cheese_0_y` | 5 × 5 µm | Square holes. |
| `delta_x`, `delta_y` | 20 µm | Centre-to-centre hole pitch. |
| `cheese.edge_nocheese` | 0.5 mm | Base renderer margin. The 10 mm wafer export enforces the 547.5 µm four-edge margin on the finished GDS. |

Layers of a fab GDS: 1 / 100 cheesed ground (cell `TOP_main_1_Cheese_100`), 1 / 10 and 1 / 11 pads,
islands, test-junction pads, CPW centre conductors. Without cheesing the ground is on 1 / 0.

## 5. Units and cell cleanup

qiskit-metal writes mm; `_prune_gds` re-reads the file with
`units='convert'` and writes it in **µm** (unit 1e-6, precision 1e-9) like the lab's files.
(gdspy's default `units='skip'` keeps the mm numbers but labels them µm: the chip then shows up
10 µm wide in KLayout.) Cells not referenced from `TOP`, and the dangling reference to
`Cheese_diff` that fabricate mode leaves, are removed.

## 6. Corner test-junction structures (lab standard)

Copied from the lab's `layout3 ftf.gds`: in each chip corner an etched pocket with two floating pads
for the fab's test junctions; the junction is not drawn.

| item | size (µm) |
|---|---|
| pocket | 800 × 300, centred at (±3950, ±4250) on a 10 × 10 mm chip |
| left pad | 225 × 129 |
| right pad | 295 × 129 |
| gap between the pads | 211 |

Same orientation in all corners (wide pad on the right). The 50 µm cheese keep-out applies around
the pocket. Shared component: `TestJunction` (`pattern design/test_junction.py`).
