# Skill: GDS export

How a mask GDS is written for any Design. Design-specific export options and checks are in
`Design<N>/notes/`.

---

## 1. Rule: the GDS is junction-free

The fab draws every Josephson junction, junction array and test junction. The mask therefore holds
**no junction cell and no junction-related metal** (leads, JJ stubs, array closing ends, SQUID
blocks), unless the Design explicitly keeps landing electrodes for the fab.

* Every component with junction-related metal has a `jj_metal` option (default `'True'` for
  simulations). The GDS is built from a design made with `jj_metal=False`.
* `export_gds` empties the `junction` qgeometry table, so no JJ library is needed.
* A new component with junction-related metal must get a `jj_metal` option too.
* Q3D / HFSS render the design directly (simulation geometry); no simulation GDS is written.

## 2. How to export

1. Before exporting, ask whether the four corner test-junction structures should be included,
   unless the user already chose for this export. Omitted: remove both the pads and the etched
   pockets. Included: export pads and pockets; the junction itself is left to the fab.
2. Use only the shared code in `pattern design/gds_tools.py`: `export_gds`, `apply_fab_gds_options`,
   `plot_gds`, `junction_window`. Never copy these into a Design.
3. `export_gds` edits the design in place (junction table emptied): build a new design for each file.
4. Write two files: without cheesing and with cheesing (`fab=True`). File names: `skill.md` section 14.

## 3. Checks after an export

1. Open the GDS in KLayout: units µm, chip outline correct (10 × 10 mm chip: −5000 … +5000 µm).
2. Junction sites are empty gaps (or hold only the landing electrodes the Design asks for).
3. Re-export after any component change; commit GDS and figures together.

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
gds_render.options.cheese.update(dict(cheese_0_x=0.005, cheese_0_y=0.005, delta_x=0.05, delta_y=0.05))
```

Design units are mm, so 0.005 = 5 µm and 0.05 = 50 µm.

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
| `delta_x`, `delta_y` | 50 µm | Hole pitch. |
| `cheese.edge_nocheese` | 0.5 mm | No holes in the outer 0.5 mm of the chip. |

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
