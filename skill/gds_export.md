# Skill: GDS export (Design1 – Design3)

How the mask GDS of the FTF designs is written. **Rule: the GDS contains no junction-related
cells.** The fab draws the FX JJs and arrays, TC Dolan junction patterns, and any test junctions.
Design3 retains both TC SQUID landing electrodes in the base-metal mask.

---

## 1. What "junction-free" means

| component (file in `SQL HFSS/pattern design/`) | left out with `jj_metal='False'` | kept |
|---|---|---|
| `MyFluxonium` (`my_fluxonium.py`, Design1 / 3) | 2 µm lead line, the metal closing end of the U (`array_end`) | `pad_top`, `pad_bot` (with the 30 µm stems), pocket, flux U-loop, readout pad, charge line |
| `MyTmonCoupler` (`my_tmon.py`, Design1 / 3) | SQUID block (island electrode, loop hole, ground electrode `squid_ground`) | bar, paddles, stem; the stem ends in its etch window |
| `FluxoniumPocket` (`fluxonium_withJJ.py`, Design2) | JJ pads `JJpad_up` / `JJpad_bottom` (arm stubs already off: `show_l_arms='False'`) | round pads, pocket, lines |
| `Tmon_c` (`tmon_c.py`, Design2) | SQUID pads `JJpad_bottom` / `JJpad_up` | stem, arms, palms, etch |
| `TestJunction` / corner test structures | junction (never drawn) | pocket + two pads |

On top of that `export_gds` empties the `junction` qgeometry table, so no junction cell is
placed (no JJ library is needed).

Default of every component is `jj_metal='True'`.
* Design1 / Design3: Q3D / HFSS keep the full geometry (`build_design()`); only the GDS is built
  with `build_design(jj_metal=False)`.
* Design2: `build_design()` itself defaults to `jj_metal=False`, so plots, GDS and Q3D all use
  the junction-free version (`build_design(jj_metal=True)` = lab pads).

Expected look (FX1 of Design1 / 3, `fig/layout/ftf_layout_fx1.png`; Design2: `ftf_layout_fxjj.png`, the two Q1 stems with an empty gap): white pocket with the
readout pad, `pad_top`, `pad_bot` and the flux U-loop; nothing between the pads. Coupler (check in KLayout): Design3 keeps both TC SQUID electrodes with two open JJ sites.
The other junction-free layouts may use an empty etched window.

## 2. How to export

Before exporting, ask whether the four corner test-junction structures should be included,
unless the user has already specified a choice for this export. With the structures omitted,
remove both the test pads and their etched pockets. With the structures included, export the
pads and pockets; leave the test junction itself for fabrication.

Design3's CLI prompts when no choice is supplied. After confirmation, pass
`--test-junctions yes` or `--test-junctions no`. Direct Python exports of Design3 should
pass `include_test_junctions=True` or `False` explicitly to `build_design`, and use
`jj_metal=False, include_tc_landing_pads=True` to retain the TC landing electrodes.
Enable `fab_contacts=True` only for the fab GDS.

```
cd "Design1/pattern design"      # or Design2 / Design3
& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I ftf_layout.py
```

In Python:

```python
from ftf_layout import build_design, export_gds      # export_gds comes from pattern design/gds_tools.py
# Design3 example after the user chooses to omit corner test structures.
options = dict(jj_metal=False, include_tc_landing_pads=True, include_test_junctions=False)
export_gds(build_design(**options), '../gds/ftf_layout.gds')
export_gds(build_design(**options, fab_contacts=True), '../gds/ftf_layout_fab.gds', fab=True)
```

Code: `SQL HFSS/pattern design/gds_tools.py` (shared by Design1 – 3): `export_gds`,
`apply_fab_gds_options`, `plot_gds`, `junction_window`. Never copy these into a
Design's `ftf_layout.py`.

`export_gds` edits the design in place (junction table emptied): build a new design for each
file.

Every Design writes `gds/ftf_layout.gds` (junction-free, no cheesing) and
`gds/ftf_layout_fab.gds` (junction-free, cheesed). The figures made from them (names, windows,
format, how to add a view) are in **`layout_figures.md`**.

`Design2/gds/layout3 ftf.gds` is the lab's reference layout (not generated). No simulation GDS
is written: Q3D / HFSS render the design directly (`build_design(for_process=True)`).

## 3. Checks after an export

1. Open the GDS in KLayout: units µm, chip −5000 … +5000 µm.
2. Check the FX sites for empty junction and array regions. On Design3, check that both
   TC landing electrodes are present, each reaches its own base conductor, and both JJ
   gaps remain open. The base-metal mask uses layer 1.
3. Rerun `python ftf_layout.py` after any component change; commit GDS + figures together.

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
| `short_segments_to_not_fillet` | False | Fillet every corner, even on short segments (all routes have segments ≥ 2 × 99 µm). |
| `tolerance` | 1e-6 m = 1 µm | Tolerance for the filleted path outlines. |
| `fabricate` | True | Write only the final mask; drop intermediate cells (keep-out, hole cell, uncheesed ground). |
| `negative_mask` | `{main: [2]}` | Layer 2 would be negative. Layer 1 (all metal here) stays **positive**. |
| `no_cheese.buffer` | 50 µm | No holes within 50 µm of any etched feature. |
| `no_cheese.view_in_file` | `{1: False, 2: False}` | Keep-out region not written to the file. |
| `cheese.view_in_file` | `{1: True, 2: False}` | Cheese the ground on layer 1 only. |
| `cheese_0_x`, `cheese_0_y` | 5 × 5 µm | Square holes. |
| `delta_x`, `delta_y` | 50 µm | Hole pitch. |
| `cheese.edge_nocheese` | 0.5 mm | No holes in the outer 0.5 mm of the chip. |

Layers of `ftf_layout_fab.gds`: 1 / 100 cheesed ground (cell `TOP_main_1_Cheese_100`),
1 / 10 and 1 / 11 pads, islands, test-junction pads, CPW centre conductors.
`ftf_layout.gds`: ground on 1 / 0.

## 5. Units and cell cleanup

qiskit-metal writes mm; `_prune_gds` re-reads the file with
`units='convert'` and writes it in **µm** (unit 1e-6, precision 1e-9) like the lab's files.
(gdspy's default `units='skip'` keeps the mm numbers but labels them µm: the chip then shows up
10 µm wide in KLayout.) Cells not referenced from `TOP`, and the dangling reference to
`Cheese_diff` that fabricate mode leaves, are removed.

## 6. Test-junction structures (chip corners, Design1 / 3)

Copied from the lab's `layout3 ftf.gds`: in each corner an etched pocket with two floating pads
for the fab's test junctions; the junction is not drawn.

| item | size (µm) |
|---|---|
| pocket | 800 × 300, centred at (±3950, ±4250) |
| left pad | 225 × 129 |
| right pad | 295 × 129 |
| gap between the pads | 211 |

Same orientation in all corners (wide pad on the right). Code: `TEST_JJ` in `ftf_layout.py`
(Design1: `Rectangle` components, Design3: `TestJunction`). The 50 µm cheese keep-out applies
around the pocket.
