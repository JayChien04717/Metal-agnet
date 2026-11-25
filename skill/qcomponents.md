# Shared QComponents (`pattern design/`)

Shared `QComponent` library used across chip layouts. Chip assemblies stay in each Design's `layout/`; Design-specific component forks stay in `Design<N>/layout/components/`. Coding conventions: [skill.md](skill.md).

## Folders and imports

| Folder | Contents |
|---|---|
| `pattern design/` | General parts: ports, launch pad, test-junction coupon, `metal_width.py` |
| `fluxonium/` | Fluxonium qubits and fab variants |
| `transmon/` | Transmon / tunable coupler, shared flux-bias line |
| `resonator/` | Resonator and Purcell-channel parts (standard channel: `PurcellReadout`) |
| `marker/` | Die / alignment markers and ID labels |
| `ebeam_junction/` | E-beam junction masks ([ebeam_junction.md](ebeam_junction.md)) |

Modules use flat names (e.g. `from my_tmon_dolan import MyTmonCouplerDolan`), so all folders above must be on `sys.path`: call `project_paths.configure(__file__)` (or `project_paths.pattern_design_paths()`). Never insert only `pattern design/` itself.

## Rules for every component

1. **ASCII sketch** of the geometry (orientation 0) with pins and option names in the module or class docstring.
2. **Minimum metal width 4 µm** ([chip_drc.md C12a](chip_drc.md)): enforce every option-defined ground strip with `metal_width.require_strip()` in `make()`. Exempt: resonator parts (including `CoupledLineTee` / coupled-line tees, readout & Purcell resonators, resonator capacitors) and the die-ID label; junction metal is outside this rule.
3. **Mask vs. simulation:** `jj_metal='False'` drops junction metal (junction-free mask); fab variants (`*_fab.py`) add contact notches only (create `*_fab.py` only in Stage 3 when preparing fab exports, not during Stage 1 draft layout); `flux_as_ground='True'` merges a shorted flux line into ground for Q3D.
4. **Labels:** die / design ID height `LABEL_SIZE_UM` = 150 µm (`marker/xmark.py`).
5. **Post-edit check:** in Stage 1, verify with `metal_width.require_strip()` / `check_component()`; in Stage 3 (fab mask export), run `python -I code/layout/check_metal_width.py <mask>.gds` (or `design_exemptions(design)` + `check()`).
6. **Keep-out polygons (markers and coupons):** marker and coupon `QComponent`s (`marker/*`, `junction_test_coupon.py`) must expose their keep-out regions (e.g. a `keepouts()` method returning Shapely polygons in design units, or named ground-opening polygons) so QGeometry-level clearance checks (`check_keepouts`, [layout_structure.md](layout_structure.md)) can verify route clearance in Stage 1 before GDS export.

## Component inventory

| File | Class / helpers | Notes |
|---|---|---|
| `pcbport.py`, `wafer_pcbport.py` | `pcb_ports`, `wafer_pcb_ports` | PCB anchors; standard 8-port wafer side ports at ±2 mm |
| `bondpad.py` | `BondPad`, `add_bond_pad` | [Standard wire-bond launch pad](gds_export.md#wire-bond-launch-pad-standard) |
| `junction_test_coupon.py` | `HookedTestJunction` | Corner / standalone test-junction coupon: floating pads, landing stems, optional Dolan recesses (F/T styles; [gds_export.md §6](gds_export.md#6-corner-test-junction-coupons)) |
| `metal_width.py` | `require_strip`, `thin_metal`, `check_component` | 4 µm minimum-metal check |
| `fluxonium/fluxonium_nthu_dolan.py` | `MyFluxoniumDolan`, `FluxoniumArrayJJDolan` | Planar (NTHU) fluxonium (charge strip `cl_sep - cpw_gap` ≥ 4 µm) and its simulation model (`flux_as_ground`, arms/closure array models) |
| `fluxonium/fluxonium_nthu_dolan_fab.py` | `MyFluxoniumDolanFab` | Adds Dolan contact notches; `contact_hooks()` returns placed notches and array direction |
| `fluxonium/fluxonium_delft_dolan.py` | `FluxoniumDelftDolan` | **Stage 1 draft (electrically unvalidated):** standalone lab (Delft / Andersen) round-pad fluxonium (`flux_line_options.signal_leg` 'top' / 'bottom'): pocket, round pads on stems, left-wall U-loop flux line, pill readout pad, offset charge line (defaults reproduce the lab-mask Q1 base metal). No Manhattan JJ drawing; JJ + array = NTHU hairpin (`MyFluxoniumDolan._loop`, towards the flux line); `jj_metal='False'` drops leads and hairpin end. Check: `code/qubit/preview_fluxonium_delft.py` |
| `fluxonium/fluxonium_delft_dolan_fab.py` | `FluxoniumDelftDolanFab` | `FluxoniumDelftDolan` with `jj_metal='False'` plus 2 × 4 µm Dolan contact notches in the facing stem ends; `contact_hooks()` as `MyFluxoniumDolanFab` |
| `ebeam_junction/fluxonium_junction_dolan.py` | `FluxoniumJunctionDolan` | E-beam mask: small JJ + JJ array + funnel contacts; finger always left; area set by `jj_cross_electrode_width` ([ebeam_junction.md](ebeam_junction.md), [doc](../pattern%20design/ebeam_junction/fluxonium_junction_dolan.md)) |
| `transmon/my_tmon_dolan.py` | `MyTmonCouplerDolan` | Tunable transmon / coupler; flux & charge ground strips ≥ 4 µm; optional readout pad via `readout_line_options.make_rol` |
| `transmon/my_tmon_dolan_fab.py` | `MyTmonCouplerDolanFab` | Adds Dolan contact notches; `contact_hooks()` returns placed island/ground notches |
| `transmon/tmon_c.py` | `Tmon_c`, `Tmon_c_B2` | Legacy transmon coupler |
| `ebeam_junction/transmon_junction_dolan.py` | `TransmonJunctionDolan` | E-beam mask: Dolan SQUID or single JJ; JJ2 translated by `squid_spacing` (no mirror); per-branch area via `jj1_/jj2_cross_electrode_width` ([doc](../pattern%20design/ebeam_junction/transmon_junction_dolan.md)) |
| `transmon/fluxline.py` | `geometry`, `B2_FLUX_OPTIONS` | L / T / U flux line; raises if ground strip to SQUID pocket < 4 µm |
| `resonator/hammer_geometry.py` | `hammer_shapes` | Hammer heads + pocket shared by `PurcellReadout` |
| `resonator/purcell_readout.py` | `PurcellReadout`, `hammer_coupon` | **Standard Purcell + readout pair:** λ/4 filter + λ/4 resonator joined only by a hammer; opposite-direction routing; `strap_points()`, `parallel_runs()`, `length_report`, `draw_mode='coupon'`, `short_ends='open'` ([skill](purcell_resonator_design.md), [method note](../note/purcell_filter_method.md#6-purcellreadout-layout-hammer-only-coupling)) |
| `resonator/enclosing_cap.py` | `EnclosingCap` | Hugger capacitor (C_in / C_p) |
| `resonator/tap_junction.py` | `TapJunction` | CPW T junction on feedline |
| `resonator/res_end.py` | `ResEnd` | Open resonator end (pin only) |
| `marker/xmark.py` | `XMark`, `add_corner_marks` | 400 µm etched window, 250 × 10 µm cross, etched label (`label_gap` 30 µm) |
| `marker/nthu_tiau_marker.py` | `TiAuChipMarker`, `TiAuWaferMarker` | NTHU TiAu e-beam markers: 4 chip sites per die (layer 2) / 8 wafer sites (layer 3), fixed photomask coordinates, optional 50 µm ground opening ([note](../note/tiau_marker_alignment.md)) |
| `marker/tsri_a12_mark.py` | `TSRIA12Mark` | TSRI A12 chip marks: 4 fine-cross groups + 4 coarse crosses within ±4.8 mm (layer 2), 50 µm ground opening per group; corner L marks are street artwork ([fab rules](fabrication_output.md#tsri-a12-marks)) |

Related shared tools (in `code/`): `code/qubit/jj_resistance_model.py`, `code/ebeam/check_ebeam_junctions.py`, `code/resonator/purcell_readout_model.py`, `code/resonator/cpw_pair_cap2d.py`, `code/layout/gds_tools.py`, `code/layout/check_metal_width.py`, `code/layout/check_keepouts.py`, `code/layout/route_length_solver.py`, and `code/fab/check_marker_references.py` (marker QComponents vs. `Fabrication output/Template/` references).
