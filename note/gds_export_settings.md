# GDS export settings (fabrication / cheesing)

Used for `gds/ftf_layout_fab.gds`. In the code they live in
`apply_fab_gds_options()` in `pattern design/ftf_layout.py`; `export_gds(..., fab=True)`
applies them.

```python
gds_render = design.renderers.gds
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

## What each setting does

Design units are mm, so 0.005 = 5 µm and 0.05 = 50 µm.

| option | value | meaning |
|---|---|---|
| `short_segments_to_not_fillet` | False | Fillet every corner, even on short segments. All routes here have segments ≥ 2 × 99 µm anyway. |
| `tolerance` | 1e-6 (m) = 1 µm | Tolerance for the rounded (filleted) path outlines. |
| `fabricate` | True | Write only the final mask; drop the intermediate cells (keep-out, single-hole cell, uncheesed ground). |
| `negative_mask` | `{main: [2]}` | Layer 2 would be written as a negative mask. Layer 1, where all metal of this design is, stays a **positive** mask. |
| `no_cheese.buffer` | 50 µm | Keep-out: no holes within 50 µm of any etched feature (pockets, CPW gaps, launch pads). |
| `no_cheese.view_in_file` | `{1: False, 2: False}` | Don't write the keep-out region into the file. |
| `cheese.view_in_file` | `{1: True, 2: False}` | Cheese the ground on layer 1 only. |
| `cheese_0_x`, `cheese_0_y` | 5 µm × 5 µm | Square holes (`cheese.shape` = 0). |
| `delta_x`, `delta_y` | 50 µm | Hole pitch. |
| `cheese.edge_nocheese` | 0.5 mm | No holes in the outer 0.5 mm of the chip. |

## Output (`gds/ftf_layout_fab.gds`)

| layer / datatype | content |
|---|---|
| 1 / 100 | Cheesed ground plane: ground with the 5 µm holes removed (cell `TOP_main_1_Cheese_100`, 32,018 polygons) |
| 1 / 10, 1 / 11 | Pads, qubit islands, CPW centre conductors (positive metal) |

- **Junctions:** none placed. Every junction site is left empty (`include_junctions=False`), see `export_gds`.
- **Cell cleanup:** after export, cells not referenced from `TOP` are pruned. So are
  references to the `Cheese_diff` cell, which qiskit deletes in fabricate mode but
  still references.

Other files in `gds/`:

| file | contents |
|---|---|
| `ftf_layout.gds` | Same layout without cheesing. Ground on 1/0. |
| `ftf_layout_forprocess.gds` | Simulation version: lumped-inductor gap instead of JJ + array. |

## Regenerate

```
cd "my design/pattern design"
python ftf_layout.py          # writes ftf_layout.gds, ftf_layout_fab.gds, ftf_layout_forprocess.gds
```
