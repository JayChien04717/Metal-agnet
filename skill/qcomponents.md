# Shared QComponents (`pattern design/`)

QComponents shared by Design1–6. Chip assemblies stay in each Design's `layout/`; Design-specific forks
(for example `Design3/layout/components/round_mask.py`) stay in that Design. Coding conventions: [skill.md](skill.md).
Main FTF layout: `Design3/layout/layout.py` `build_design()` (see [Design3](../Design3/README.md)); `tc_readout=True` adds the transmon readout pad and resonator, both variants use the shared `transmon/my_tmon.py`.

## Folders and imports

| Folder | Contents |
|---|---|
| `pattern design/` | General parts: ports, launch pad, test junctions, `metal_width.py` |
| `fluxonium/` | Fluxonium qubits, fab variants, junction mask drawing |
| `transmon/` | Transmon / tunable coupler, shared flux-bias line |
| `resonator/` | Resonator and Purcell-channel parts |
| `marker/` | Die / alignment markers and ID labels |

Modules keep flat names (`from my_tmon import MyTmonCoupler`), so all five folders must be on `sys.path`:
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
| `test_junction.py` | `TestJunction`, `HookedTestJunction` | Corner test pads; aligned or offset landing recesses |
| `metal_width.py` | `require_strip`, `thin_metal`, `check_component` | 4 µm rule |
| `fluxonium/my_fluxonium.py` | `MyFluxonium` | Flux line 80 / 11 µm (Design3; Design1 50 / 15 µm); charge strip `cl_sep - cpw_gap` ≥ 4 µm |
| `fluxonium/my_fluxonium_fab.py` | `MyFluxoniumFab` | Contact notches in the stem ends |
| `fluxonium/fluxonium.py` | `FluxoniumArrayJJ` | Simulation model, same loop as `MyFluxonium._loop` (JJ in the centre lead, U sides lumped L, metal U end) |
| `fluxonium/fluxonium_withJJ.py` | `FluxoniumPocket` | Design2; `u_array='True'` for the same loop |
| `fluxonium/fluxonium_junction.py` | `FluxoniumJunction` | Multilayer junction mask draft ([source and usage](../pattern%20design/fluxonium/fluxonium_junction.md)) |
| `transmon/my_tmon.py` | `MyTmonCoupler` | Flux and charge ground strips ≥ 4 µm; readout pad only with `readout_line_options.make_rol` |
| `transmon/my_tmon_fab.py` | `MyTmonCouplerFab` | Dolan contact notches |
| `transmon/tmon_c.py` | `Tmon_c`, `Tmon_c_B2` | Legacy (Design2) |
| `transmon/fluxline.py` | `geometry`, `B2_FLUX_OPTIONS` | L / T / U line; B2: bar 45 µm, CPW 10/6 µm, `fbl_sep` 10.5 µm (9 µm before 2026-11-18), x +12 µm; raises if the strip to the SQUID pocket < 4 µm |
| `resonator/hammer_cap.py` | `HammerCap` | C_j |
| `resonator/enclosing_cap.py` | `EnclosingCap` | Hugger, C_in / C_p |
| `resonator/interdigital_cap.py` | `InterdigitalCap` | Finger coupler |
| `resonator/tap_junction.py` | `TapJunction` | CPW T junction on the feedline |
| `resonator/res_end.py` | `ResEnd` | Open resonator end (pin only) |
| `marker/xmark.py` | `XMark`, `add_corner_marks` | 400 µm etched window, 250 × 10 µm cross, etched label right of / below the cross, `label_gap` 30 µm |

Related tools (not QComponents): `code/general/gds_tools.py` (`export_gds`, `plot_gds`, `junction_window`;
[gds_export.md](gds_export.md), [layout_figures.md](layout_figures.md)) and `code/general/check_metal_width.py`.
