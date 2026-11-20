# Shared QComponents (`pattern design/`)

QComponents shared by Design1–6. Chip assemblies stay in each Design's `layout/`; Design-specific forks
(for example `Design3/layout/components/round_mask.py`) stay in that Design. Coding conventions: [skill.md](skill.md).
Main FTF layout: `Design3/layout/layout.py` `build_design()` (see [Design3](../Design3/README.md)); `tc_readout=True` adds the transmon readout pad and resonator, both variants use the shared `transmon/my_tmon.py`.

## Folders and imports

| Folder | Contents |
|---|---|
| `pattern design/` | General parts: ports, launch pad, test-junction coupon, `metal_width.py` |
| `fluxonium/` | Fluxonium qubits and fab variants |
| `transmon/` | Transmon / tunable coupler, shared flux-bias line |
| `resonator/` | Resonator and Purcell-channel parts (new channels: `PurcellReadout`) |
| `archive/` | Not recommended, kept importable for old Design scripts ([why](../pattern%20design/archive/README.md)) |
| `marker/` | Die / alignment markers and ID labels |
| `ebeam_junction/` | E-beam junction masks ([skill](ebeam_junction.md)) |

Modules keep flat names (`from my_tmon import MyTmonCoupler`), so all these folders must be on `sys.path`:
`project_paths.configure(__file__)` adds them, and `project_paths.pattern_design_paths()` returns them.
Never insert only `pattern design/` itself.

## Rules for every component

1. **ASCII sketch** of the geometry (orientation 0) with pins and option names in the module or class docstring.
2. **Minimum metal width 4 µm** ([chip_drc.md C12a](chip_drc.md)). Every ground strip a component's options define
   is checked with `metal_width.require_strip()` in `make()`, so an invalid option set raises instead of drawing
   a thin strip. Exempt: resonator parts and the die-ID label. Junction metal is outside this rule.
3. **Mask vs. simulation**: `jj_metal='False'` drops junction metal (junction-free mask); fab variants
   (`*_fab.py`) add contact notches only; `flux_as_ground='True'` merges a shorted flux line into ground for Q3D.
4. **Labels**: die / design ID height `LABEL_SIZE_UM` = 150 µm (`marker/xmark.py`; 100 µm before 2026-11-18).
5. After a geometry change, check the exported mask:
   `python -I code/general/check_metal_width.py <mask>.gds` (or `design_exemptions(design)` + `check()`).

## Components

| File | Class | Notes |
|---|---|---|
| `pcbport.py`, `wafer_pcbport.py` | `pcb_ports`, `wafer_pcb_ports` | PCB anchors; wafer side ports at ±2 mm |
| `bondpad.py` | `BondPad`, `add_bond_pad` | [Design6 launch standard](launch_pad.md) |
| `junction_test_coupon.py` | `HookedTestJunction` | Corner / standalone test-junction coupon: two floating pads, landing stems, optional Dolan recesses (F: facing, T: offset); placed by `DesignTestJunction/layout/build.py` `add_test_junction` ([gds_export.md §6](gds_export.md#6-corner-test-junction-coupons)) |
| `metal_width.py` | `require_strip`, `thin_metal`, `check_component` | 4 µm rule |
| `fluxonium/my_fluxonium.py` | `MyFluxonium` | Flux line 80 / 11 µm (Design3; Design1 50 / 15 µm); charge strip `cl_sep - cpw_gap` ≥ 4 µm |
| `fluxonium/my_fluxonium_fab.py` | `MyFluxoniumFab` | Contact notches in the stem ends |
| `fluxonium/fluxonium.py` | `FluxoniumArrayJJ` | Simulation model, same loop as `MyFluxonium._loop` (JJ in the centre lead, U sides lumped L, metal U end) |
| `fluxonium/fluxonium_withJJ.py` | `FluxoniumPocket` | Design2; `u_array='True'` for the same loop |
| `ebeam_junction/fluxonium_junction.py` | `FluxoniumJunction` | E-beam mask: small JJ + JJ array + funnel contacts; finger always left; `side='left'/'right'`; area set by `jj_cross_electrode_width` ([skill](ebeam_junction.md), [doc](../pattern%20design/ebeam_junction/fluxonium_junction.md)) |
| `transmon/my_tmon.py` | `MyTmonCoupler` | Flux and charge ground strips ≥ 4 µm; readout pad only with `readout_line_options.make_rol` |
| `transmon/my_tmon_fab.py` | `MyTmonCouplerFab` | Dolan contact notches |
| `transmon/tmon_c.py` | `Tmon_c`, `Tmon_c_B2` | Legacy (Design2) |
| `ebeam_junction/transmon_junction.py` | `TransmonJunction` | E-beam mask: Dolan SQUID or single JJ; JJ2 = JJ1 translated by `squid_spacing` (no mirror); per-branch area via `jj1_/jj2_cross_electrode_width` ([doc](../pattern%20design/ebeam_junction/transmon_junction.md)) |
| `archive/*_junction_layout3.py` | `FluxoniumJunction`, `TransmonJunction` (old) | **Archived** layout3 drawings (mirrored SQUID, `junction_gap`), kept for Design3 r11 / DesignTestJunction r3 ([archive README](../pattern%20design/archive/README.md)) |
| `transmon/fluxline.py` | `geometry`, `B2_FLUX_OPTIONS` | L / T / U line; B2: bar 45 µm, CPW 10/6 µm, `fbl_sep` 10.5 µm (9 µm before 2026-11-18), x +12 µm; raises if the strip to the SQUID pocket < 4 µm |
| `resonator/hammer_geometry.py` | `hammer_shapes` | Hammer heads + pocket, shared by `PurcellReadout` and the archived `HammerCap` |
| `resonator/purcell_readout.py` | `PurcellReadout`, `hammer_coupon` | **Standard Purcell + readout pair (2026-11-19).** λ/4 Purcell filter + λ/4 readout resonator joined only by a hammer; pins `cp` (C_p side, +y) and `qubit` (−y); options `mirror`, `*_first_run` ('along' / 'back'), `resonator_open_extension` (route into a side claw).
 Used by Design5 (`chip6_purcell_resonator_pair.py`), Design6 (`readout='purcell_resonator_pair'`), Design8. Lines leave the T region in opposite directions (only 2 × (`tee_half` + `lead`) side by side at `line_spacing` ≥ 400 µm). `strap_points()` gives the hammer/T air-bridge positions, `parallel_runs()` the side-by-side runs, `length_report` the filleted lengths, `draw_mode='coupon'` the Q3D hammer coupon, `short_ends='open'` the Q3D rendering. Test and validation: [Design8](../Design8/README.md); circuit model: `code/resonator/purcell_readout_model.py` (below); layout strategy: [method note section 6](../note/purcell_filter_method.md#6-purcellreadout-layout-hammer-only-coupling) |
| `resonator/enclosing_cap.py` | `EnclosingCap` | Hugger, C_in / C_p |
| `resonator/tap_junction.py` | `TapJunction` | CPW T junction on the feedline (n0 tee for C_p) |
| `archive/hammer_cap.py` | `HammerCap` | **Archived, not recommended** (old separately routed pair; the capacitance matrix needs extra corrections, see the archive README) |
| `archive/interdigital_cap.py` | `InterdigitalCap` | **Archived, not recommended** (finger C_j / C_in of the transmission filter) |
| `resonator/res_end.py` | `ResEnd` | Open resonator end (pin only) |
| `marker/xmark.py` | `XMark`, `add_corner_marks` | 400 µm etched window, 250 × 10 µm cross, etched label right of / below the cross, `label_gap` 30 µm |

Related tools (not QComponents): `code/qubit/jj_resistance_model.py` (EJ ↔ Rn ↔ actual area ↔ e-beam mask widths, `size_from_ej`, `squid_options`, `squid_sizes` (from EJ_sum and asymmetry), `ej_from_size`, `finger_from_ej`; [skill](ebeam_junction.md); checks `code/general/check_ebeam_junctions.py`), `code/resonator/purcell_readout_model.py` (Purcell readout circuit model: length ↔ frequency `tapped_lambda4_f`, `length_for_f`, `lambda4_f`, `lambda2_f`; hammer C_j fit to Q3D `hammer_cj`, `gap_for_cj`, `height_for_cj`; J / κ `j_hammer`, `kappa_hanger`, `kappa_readout`, `coupling`; parallel-run parasitic J `parallel_j` with the 2-D table `CPW_PAIR_2D`), `code/resonator/cpw_pair_cap2d.py` (2-D CPW pair capacitance), `code/general/gds_tools.py` (`export_gds`, `plot_gds`, `junction_window`;
[gds_export.md](gds_export.md), [layout_figures.md](layout_figures.md)) and `code/general/check_metal_width.py`.
