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

## Units

All GDS files in `gds/` are written in **µm** (GDS unit 1e-6, precision 1e-9); the chip spans
−5000 … +5000 µm. qiskit-metal writes mm; `_prune_gds` re-reads the file with
`units='convert'` and writes it in µm. (Before 2026-10-21 the re-read used gdspy's default
`units='skip'`, which kept the mm numbers but labelled them µm, so the chip showed up as 10 µm
wide in KLayout.)

## Test-junction structures (chip corners)

Copied from the lab's `layout3 ftf.gds`: in each chip corner an etched pocket with two
floating pads, for the fab's test junctions. The junction between the pads is **not** drawn
(same rule as the qubit junctions); the fab adds it.

| item | size (µm) |
|---|---|
| pocket (no ground) | 800 × 300, centred at (±3950, ±4250) |
| left pad | 225 × 129 |
| right pad | 295 × 129 |
| gap between the pads | 211 |

Same orientation in all four corners (wide pad on the right). In the code: `TEST_JJ` in
`ftf_layout.py`, drawn as `Rectangle` components `testjj_<LL|LR|UL|UR>_pocket` (subtract)
and `_padL`, `_padR` (metal). The 50 µm cheese keep-out applies around the pocket.

## Output (`gds/ftf_layout_fab.gds`)

| layer / datatype | content |
|---|---|
| 1 / 100 | Cheesed ground plane: ground with the 5 µm holes removed (cell `TOP_main_1_Cheese_100`, 32,758 polygons) |
| 1 / 10, 1 / 11 | Pads, qubit islands, test-junction pads, CPW centre conductors (positive metal) |

- **Junctions:** none placed. Every junction site is left empty (`include_junctions=False`), see `export_gds`.
- **Cell cleanup:** after export, cells not referenced from `TOP` are pruned. So are
  references to the `Cheese_diff` cell, which qiskit deletes in fabricate mode but
  still references.

Other files in `gds/`:

| file | contents |
|---|---|
| `ftf_layout.gds` | Same layout without cheesing (with the test-junction pockets). Ground on 1/0. |
| `ftf_layout_forprocess.gds` | Simulation version: lumped-inductor gap instead of JJ + array. |

Plots of the fab file: `fig/fab_cheese_chip.png`, `fab_cheese_testjj.png`, `fab_cheese_holes.png`.

## Regenerate

```
cd "my design/pattern design"
python ftf_layout.py          # writes ftf_layout.gds, ftf_layout_fab.gds, ftf_layout_forprocess.gds
```
