# Shared QComponents (`pattern design/`)

Shared `QComponent` library used across chip layouts. Chip assemblies stay in each Design's `layout/`; Design-specific component forks stay in `Design<N>/layout/components/`. Coding conventions: [skill.md](skill.md).

## Folders and imports

| Folder | Contents |
|---|---|
| `pattern design/` | General parts: ports, launch pad, test-junction coupon, `metal_width.py` |
| `fluxonium/` | Fluxonium qubits and fab variants |
| `transmon/` | Transmon / tunable coupler, shared flux-bias line |
| `resonator/` | Resonator and Purcell-channel parts (standard channel: `PurcellReadout`) |
| `archive/` | Deprecated parts kept importable for legacy scripts ([why](../pattern%20design/archive/README.md)) |
| `marker/` | Die / alignment markers and ID labels |
| `ebeam_junction/` | E-beam junction masks ([ebeam_junction.md](ebeam_junction.md)) |

Modules use flat names (e.g. `from my_tmon import MyTmonCoupler`), so all folders above must be on `sys.path`: call `project_paths.configure(__file__)` (or `project_paths.pattern_design_paths()`). Never insert only `pattern design/` itself.

## Rules for every component

1. **ASCII sketch** of the geometry (orientation 0) with pins and option names in the module or class docstring.
2. **Minimum metal width 4 µm** ([chip_drc.md C12a](chip_drc.md)): enforce every option-defined ground strip with `metal_width.require_strip()` in `make()`. Exempt: resonator parts and the die-ID label; junction metal is outside this rule.
3. **Mask vs. simulation:** `jj_metal='False'` drops junction metal (junction-free mask); fab variants (`*_fab.py`) add contact notches only (create `*_fab.py` only in Stage 3 when preparing fab exports, not during Stage 1 draft layout); `flux_as_ground='True'` merges a shorted flux line into ground for Q3D.
4. **Labels:** die / design ID height `LABEL_SIZE_UM` = 150 µm (`marker/xmark.py`).
5. **Post-edit check:** in Stage 1, verify with `metal_width.require_strip()` / `check_component()`; in Stage 3 (fab mask export), run `python -I code/general/check_metal_width.py <mask>.gds` (or `design_exemptions(design)` + `check()`).

## Component inventory

| File | Class / helpers | Notes |
|---|---|---|
| `pcbport.py`, `wafer_pcbport.py` | `pcb_ports`, `wafer_pcb_ports` | PCB anchors; standard 8-port wafer side ports at ±2 mm |
| `bondpad.py` | `BondPad`, `add_bond_pad` | [Standard wire-bond launch pad](gds_export.md#wire-bond-launch-pad-standard) |
| `junction_test_coupon.py` | `HookedTestJunction` | Corner / standalone test-junction coupon: floating pads, landing stems, optional Dolan recesses (F/T styles; [gds_export.md §6](gds_export.md#6-corner-test-junction-coupons)) |
| `metal_width.py` | `require_strip`, `thin_metal`, `check_component` | 4 µm minimum-metal check |
| `fluxonium/my_fluxonium.py` | `MyFluxonium` | Planar fluxonium; charge strip `cl_sep - cpw_gap` ≥ 4 µm |
| `fluxonium/my_fluxonium_fab.py` | `MyFluxoniumFab` | Adds Dolan contact notches; `contact_hooks()` returns placed notches and array direction |
| `fluxonium/fluxonium.py` | `FluxoniumArrayJJ` | Simulation model matching `MyFluxonium._loop` (centre JJ + two lumped-L U sides + metal U closing end) |
| `fluxonium/fluxonium_withJJ.py` | `FluxoniumPocket` | Round-pad fluxonium variant; `u_array='True'` enables the U-array loop |
| `ebeam_junction/fluxonium_junction.py` | `FluxoniumJunction` | E-beam mask: small JJ + JJ array + funnel contacts; finger always left; area set by `jj_cross_electrode_width` ([ebeam_junction.md](ebeam_junction.md), [doc](../pattern%20design/ebeam_junction/fluxonium_junction.md)) |
| `transmon/my_tmon.py` | `MyTmonCoupler` | Tunable transmon / coupler; flux & charge ground strips ≥ 4 µm; optional readout pad via `readout_line_options.make_rol` |
| `transmon/my_tmon_fab.py` | `MyTmonCouplerFab` | Adds Dolan contact notches; `contact_hooks()` returns placed island/ground notches |
| `transmon/tmon_c.py` | `Tmon_c`, `Tmon_c_B2` | Legacy transmon coupler |
| `ebeam_junction/transmon_junction.py` | `TransmonJunction` | E-beam mask: Dolan SQUID or single JJ; JJ2 translated by `squid_spacing` (no mirror); per-branch area via `jj1_/jj2_cross_electrode_width` ([doc](../pattern%20design/ebeam_junction/transmon_junction.md)) |
| `transmon/fluxline.py` | `geometry`, `B2_FLUX_OPTIONS` | L / T / U flux line; raises if ground strip to SQUID pocket < 4 µm |
| `resonator/hammer_geometry.py` | `hammer_shapes` | Hammer heads + pocket shared by `PurcellReadout` |
| `resonator/purcell_readout.py` | `PurcellReadout`, `hammer_coupon` | **Standard Purcell + readout pair:** λ/4 filter + λ/4 resonator joined only by a hammer; opposite-direction routing; `strap_points()`, `parallel_runs()`, `length_report`, `draw_mode='coupon'`, `short_ends='open'` ([skill](purcell_resonator_design.md), [method note](../note/purcell_filter_method.md#6-purcellreadout-layout-hammer-only-coupling)) |
| `resonator/enclosing_cap.py` | `EnclosingCap` | Hugger capacitor (C_in / C_p) |
| `resonator/tap_junction.py` | `TapJunction` | CPW T junction on feedline |
| `resonator/res_end.py` | `ResEnd` | Open resonator end (pin only) |
| `marker/xmark.py` | `XMark`, `add_corner_marks` | 400 µm etched window, 250 × 10 µm cross, etched label (`label_gap` 30 µm) |
| `archive/*` | `HammerCap`, `InterdigitalCap`, `*_junction_layout3.py` | **Archived, not recommended** for new designs ([archive README](../pattern%20design/archive/README.md)) |

Related shared tools (in `code/`): `code/qubit/jj_resistance_model.py`, `code/general/check_ebeam_junctions.py`, `code/resonator/purcell_readout_model.py`, `code/resonator/cpw_pair_cap2d.py`, `code/general/gds_tools.py`, and `code/general/check_metal_width.py`.
