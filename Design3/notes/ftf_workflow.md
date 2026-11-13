# FTF workflow (Design1 – Design3): layout, GDS export, figures

Design-specific part of the FTF chips (fluxonium – tunable transmon coupler – fluxonium), shared by
Design1, Design2 and Design3. The general rules are in `../../skill/` (`skill.md`, `gds_export.md`,
`layout_figures.md`); this file only holds what is particular to `ftf_layout.py` and its components.

## 1. Layout script `ftf_layout.py`

The chip is built by `pattern design/ftf_layout.py` (no notebook). `build_design()`
returns the design; `export_gds()` writes `gds/`; `plot_gds()` draws figures from the
exported GDS (both in the shared `pattern design/gds_tools.py`, see `gds_export.md`), so fillets and junction gaps look exactly as they will on the mask.

### 1.1 Where the knobs are
* `FX_OPTIONS` (shared fluxonium geometry) and `FX_PER_QUBIT` (per-qubit overrides, e.g. the readout gap that sets g for each resonator frequency).
* `TMON_OPTIONS`: coupler. `pad_arm_length` and the coupler height are **derived** from the FX pads (`PADDLE_SEP` = paddle-to-pad gap), so changing the FX pad size moves the coupler with it.
* `RESONATORS`: per resonator `f`, `kappa`, `C_pad` (Q3D readout-pad capacitance), `C_extra` (full-circuit Q3D coupling not in the tee fit) and `hfss_scale` (f_HFSS / f_target). The tee coupling length and the total length are computed from these.
* `FX_X`, `Y0`, `FEED_Y`: qubit positions and feedline height.
* `build_design(for_process=True, fx_options=..., tmon_options=..., fx_x=..., fx1_options=..., fx2_options=...)` overrides any of these for tuning sweeps.

### 1.2 Component conventions used here
* Draw in local coordinates, then mirror / rotate / translate with one helper (`_place`), and create pins through a helper that swaps the point order when mirrored, so normals stay outward.
* Check every pin's `middle` and `normal` after a change (`design.components[name].pins`).
* Boolean results (`draw.subtract`) can leave zero-area slivers; drop them before `add_qgeometry` (`_clean` in `my_tmon.py`), otherwise Ansys fails with *Cannot create rectangle with zero length/width*.
* Junction-table lines: draw all junctions of a device in the same direction (Dolan bridges need parallel junctions).
* **One fluxonium simulation model for every run** (Q3D, HFSS eigenmode, mutual inductance): `FluxoniumArrayJJ` (`fluxonium.py`). Same pads / stems / pocket / lines as the mask; the loop is `MyFluxonium._loop`: lead line from the stem ends with the JJ as its short middle gap at the qubit centre (`L_j`, or the port in the mutual run); U = two junction-table sides `array_top` / `array_bot` (`l_inductance / 2` each) joined by a metal closing end `array_end` (an HFSS lumped sheet is straight and must end on metal; the closing end is a floating net in Q3D, Kron-reduced by `cap_matrix.py`). No separate JJ branch. Do not make a fluxonium simulation version with a single straight inductor. Design2 (`FluxoniumPocket`) uses the same model through `u_array='True'` (`build_design(for_process=True)`). Flux lines shorted to ground are drawn as ground for Q3D (`flux_as_ground`), to avoid net-overlap errors.

### 1.3 Routing rules
* `RouteAnchors` with explicit anchor points for control lines (Pathfinder fails easily in a crowded chip); `RouteStraight` between neighbouring tees.
* Readout resonators: `RouteMixed` with one meandered segment below the tee and straight segments for the jog to the qubit; its `total_length` is not exact and fillets shorten the line, so iterate `total_length` on the **filleted** centre-line length (`filleted_length`, `path_length`) until it matches the target.
* 99 µm fillets everywhere; after building, assert every inner path segment is ≥ 2 × fillet.
* Ports: P7/P8 feedline, P1/P6 FX flux, P2/P5 FX charge, P3 TC flux, P4 TC charge; check visually that no routes cross.

### 1.4 GDS export (full rules: `gds_export.md`)
* The GDS is **junction-free**: `export_gds(build_design(jj_metal=False), file)`. `jj_metal='False'`
  drops leads / JJ stubs / array end (FX) and the SQUID block (coupler); `export_gds` empties the
  junction table. Simulations keep the default `jj_metal='True'`.
* A new component with junction-related metal must get a `jj_metal` option too.
* `fab=True` applies the lab cheesing settings; unreferenced cells are pruned (`_prune_gds`).

### 1.5 Verification after a change
* `python ftf_layout.py` runs without warnings about short segments.
* Pins / normals printed and sane; no crossing routes in `fig/layout/ftf_layout_chip.png`.
* Resonator lengths: target vs filleted length printed equal.
* Re-run the relevant Q3D / HFSS check (see `cap_matrix_analysis.md`) before quoting numbers.

## 2. GDS export of the FTF chips

### 2.1 What "junction-free" means

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

### 2.2 How to export

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

### 2.3 Checks after an export

1. Open the GDS in KLayout: units µm, chip −5000 … +5000 µm.
2. Check the FX sites for empty junction and array regions. On Design3, check that both
   TC landing electrodes are present, each reaches its own base conductor, and both JJ
   gaps remain open. The base-metal mask uses layer 1.
3. Rerun `python ftf_layout.py` after any component change; commit GDS + figures together.

### 2.4 Corner test junctions in the FTF chips
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

## 3. Layout figures of the FTF chips

### 3.1 Standard figure set (`fig/layout/`, four files)

Run once per Design (Python of the `qmetal` env, see `skill.md` §9):

```
cd "D:\...\SQL HFSS\Design1\pattern design"        # Design2, Design3 the same
& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I ftf_layout.py
```

This writes `gds/ftf_layout.gds`, `gds/ftf_layout_fab.gds` and the four figures below. Windows are
`(x_min, x_max, y_min, y_max)` in mm (Design1 / 3: `Y0 = -0.8`).

| file (`fig/layout/`) | shows | Design1 / Design3 window | Design2 window |
|---|---|---|---|
| `ftf_layout_chip.png` | whole chip | (−5, 5, −5, 5) + port labels | (−5, 5, −5, 5) |
| `ftf_layout_zoom.png` | the FTF set | (−1.0, 1.0, Y0 − 0.55, Y0 + 0.55) | (−3.2, −1.5, −2.9, −0.7) |
| `ftf_layout_fx1.png` | FX1 / Q1 with flux, charge, readout lines | (−0.95, −0.3, Y0 − 0.3, Y0 + 0.3) | (−2.85, −1.75, −2.82, −1.72) |
| `ftf_layout_fxjj.png` | FX1 / Q1 junction site (must be empty) | (−0.48, −0.32, Y0 − 0.045, Y0 + 0.045), dpi 140 | `junction_window(d, 'Q1', half=0.12, half_y=0.045)`, dpi 140 |

`junction_window(design, comp)` centres the view on the junction rows of a component. Call it
**before** `export_gds`, which empties the junction table.

No other layout pictures are made (no SQUID, die-ID or cheesed-fab plots); check those in
KLayout if needed, without saving pictures.

### 3.2 Checklist before committing FTF figures

- [ ] Made by `ftf_layout.py` from the current GDS (time stamps of GDS and PNG agree).
- [ ] Junction site empty in `ftf_layout_fxjj.png`; no red in any plot.
- [ ] `Design1/fig/layout`, `Design2/fig/layout`, `Design3/fig/layout` hold exactly the four files.
- [ ] README image links point to existing files (`Design1/3`: `ftf_layout_chip.png`, `Design2`: `ftf_layout_zoom.png`).

## 6. Other figure folders (not layout)

| folder | made by | note |
|---|---|---|
| `fig/mutual/` | `simulation/plot_mutual_z.py` (impedance); `fluxonium_arrayjj*.png` from the `__main__` of the shared `pattern design/fluxonium.py` (`plot_design`, qgeometry tables) | simulation geometry (lumped-L array), not the mask |
| `fig/zz/` | `simulation/zz_analysis.py` | ZZ maps, `../note/zz_method.md` |
