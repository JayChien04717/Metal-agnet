# Skill: Single-die GDS export, launch pads, cheesing, and test coupons

How a single-die mask GDS is built and exported for any Design. Design-specific export options, port assignments, and F/T coupon assignments belong in each Design's `README.md` / `design.json` / `notes/`.

For 4-inch wafer assembly (streets, placement IDs, tone inversion, e-beam site markers, and delivery packaging), see [fabrication_output.md](fabrication_output.md).

---

## 1. Electrical geometry and fabrication

The fab draws every Josephson junction, junction array, and test junction. The base-metal mask therefore holds **no junction cell and no junction-related metal** (leads, JJ stubs, array closing ends, SQUID blocks), unless a Design explicitly keeps landing electrodes for the fab. Drawn e-beam junction layers (11–14, [ebeam_junction.md](ebeam_junction.md)) are added only when chosen in step 3 of Section 2.1.

* Every component with junction-related metal has a `jj_metal` option (default `'True'` for simulations). The mask GDS is built with `jj_metal=False`.
* `export_gds` empties the `junction` QGeometry table, so no JJ library is needed.
* Q3D and HFSS render the unnotched simulation `QComponent` directly from QGeometry (the junction table provides lumped elements). A GDS export cannot encode those boundary conditions—never import a GDS into a solver.
* A companion simulation/reference GDS (`*_sim.gds`) preserves solver geometry on layer 1 and annotates lumped junction sheets on layer 50/0 (record and figure source only).
* Fabrication `QComponent` variants (`*_fab.py`, enabled via `fab_contacts=True`) add Dolan contact hooks (2 µm wide × 4 µm deep; transmon ground-side cuts use a signed +2 µm local x offset). Confirm clear junction gaps and that hooks do not sever the landing pads. Resonator-only chips need no junction hooks.

---

## 2. How to export

There is one geometry source: the Design's `layout/` `build_design()`. Simulations render it directly into Q3D / HFSS (`jj_metal=True`, lumped junctions, no hooks, no cheese). GDS files are a separate render of the **same** `build_design()`. Never edit a GDS by hand.

> **Stage gate reminder ([skill.md §2.0](skill.md#20-three-stage-gate-rule-do-not-jump-ahead)):** During **Stage 1 (draft layout & visual review)**, export **only** the clean reference `*_sim.gds` so `code/layout/plot_design_figs.py` can render `fig/chip.png` for user inspection. Do **not** build `*_fab.py` hooks, place e-beam junctions, export `*_fab.gds`, or run fab mask checks until **Stage 3** (after the user has reviewed the layout figures, confirmed simulation parameters, and requested fabrication output).

### 2.1 Fab GDS workflow (Stage 3 only)

<a id="21-fab-gds-workflow"></a>

Produce a fab GDS in this order—do not skip or reorder steps:

1. **Simulation parameters confirmed:** the Design's `design.json` records adopted runs and `final_selection` ([parameter_provenance.md](parameter_provenance.md)). Use those exact `build_design()` options with no geometry retuning. If unconfirmed, stop and ask.
2. **Swap to fab components:** build with `jj_metal=False` and `fab_contacts=True` (plus any Design-specified landing-pad option). Only this swap differs from the simulation geometry.
3. **Ask the user** about each item below (unless already specified for this export) and record the answers:

   | Item | Question | Rules |
   |---|---|---|
   | Single die or wafer | One 10 mm die, or a 4-inch wafer assembly (which grid)? | Section 3, [fabrication_output.md](fabrication_output.md) |
   | E-beam site | NTHU, TSRI, or NCKU? | [fabrication_output.md §2](fabrication_output.md#2-alignment-markers-by-e-beam-site) |
   | Alignment markers | Site markers and bare rectangles (NTHU TiAu wafer + chip; TSRI pending; NCKU none) | [fabrication_output.md §2](fabrication_output.md#2-alignment-markers-by-e-beam-site) |
   | E-beam junction layers | Draw junctions (`FluxoniumJunctionDolan` / `TransmonJunctionDolan`) or leave empty for fab? If drawn: per-JJ EJ and layer map | [ebeam_junction.md](ebeam_junction.md) |
   | Corner test junctions | Include the four F/T coupons? Which corner gets F / T? | Section 6 |
   | Die ID / design ID | Label text (must not start with a component name) | Section 3 |

4. **Render** with `export_gds(design, path, fab=True)` from `code/layout/gds_tools.py`:
   - **Cheese goes in last**, after hooks, coupons, junction layers, and marker openings are final, so its 50 µm keep-out sees every etch.
   - Single die: `fab=True` cheeses the ground (Section 4).
   - Wafer: render with `cheese=False` and rebuild the ground grid on the assembled wafer (`wafer_cheese_grid.replace_wafer_cheese`, Section 3).
5. **Check** per Section 3, then package deliveries per [fabrication_output.md](fabrication_output.md).

**E-beam junction layers in the render:** when drawn junctions are chosen in step 3, place them inside `build_design()` via [`code/ebeam/ebeam_junction_placement.py`](../code/ebeam/ebeam_junction_placement.py) `add_ebeam_junctions()`, reading placed hooks from `MyFluxoniumDolanFab.contact_hooks()` / `MyTmonCouplerDolanFab.contact_hooks()` on layers 11–14. Render with `export_gds(design, path, fab=True, datatype0_layers=(11, 12, 13, 14))` so positive shapes write as datatype 0. A simulation build must refuse `ebeam_junctions=True`.

### 2.2 Code rules

1. Use only the shared functions in `code/layout/gds_tools.py` (`export_gds`, `apply_fab_gds_options`, `plot_gds`, `junction_window`). Never copy them into a Design.
2. `export_gds` mutates the design in place (emptying the junction table), so build a fresh design for each output file.
3. In **Stage 1 (draft layout)**, write only the unnotched reference / simulation GDS (`*_sim.gds`) required by `plot_design_figs.py`. In **Stage 3 (fab export)**, write the matching pair: `*_sim.gds` and `*_fab.gds`.

### 2.3 Which GDS files may be written

Write a GDS only when something downstream consumes it. Every GDS left on disk must be one of these:

| Role | Where | When | Consumer |
|---|---|---|---|
| `<Design>_<date>[_<variant>]_rN_sim.gds` | `Design<N>/gds/` | Stage 1 onward | layout figures, record of solver geometry |
| `<Design>_<date>[_<variant>]_rN_fab.gds` | `Design<N>/gds/` | Stage 3 only | fab / wafer assembly input |
| Delivery files (`_positive`, `_etch`, `_ebeam_<SITE>`) | `Fabrication output/<date>/` | Stage 3, on request | the fab ([fabrication_output.md §4](fabrication_output.md#4-folder-layout-and-marker-references)) |
| Fixed marker references | `Fabrication output/Template/` | never rewritten | scripts read them in place |

Rules:

1. **One file per role per revision.** `_<variant>` is allowed only for a genuinely different chip that will itself be simulated or fabricated, not for a view of the same chip.
2. **No inspection GDS.** Previews, marker-only layers, overlays, before/after diffs, single-layer extracts, and intermediate booleans are not saved as GDS. Show them as PNG figures or record the check in JSON. If a script needs an intermediate GDS, stage it in a temporary folder and delete it when the script ends (as `export_wafer_test_junctions.py` and `invert_wafer_tone.py` do). Script options that write extra GDS (e.g. `tiau_marker_overlay.py --marker-gds`) stay off unless the user asks.
3. **Write the final name directly.** Never copy or rename a GDS to make another file; two files with the same SHA256 are an error.
4. **Do not copy a GDS into another folder** to make it available; point the consuming script at the original path.
5. **Superseded revisions:** a new revision gets a new `rN`; move superseded files out of the active `gds/` folder (`gds/previous/` or delete; git keeps history) unless `design.json` still adopts them.
6. **Before ending the turn,** list every GDS written in the turn with its role, and delete anything that has no row in the table above.

Older Designs contain legacy roles (`_reference`, `_mask`, `_ebj_preview`, `_ebeamjj`). Keep the existing files, but do not create new files with these suffixes: a reference is `_sim`, a drawn e-beam junction build is `_fab`, and a preview is a PNG.

---

## 3. Checks after an export

1. **KLayout check:** units µm, chip outline exact (10 × 10 mm die: −5000 … +5000 µm).
2. **Junction sites:** empty gaps (or only explicitly requested landing electrodes).
3. **Sync:** re-export after any geometry change; commit GDS and figures together ([layout_figures.md](layout_figures.md)).
4. **Etch survival probe:** verify every etched feature survived export (e.g. `gdstk` `polygon.contain` inside each flux-line / CPW gap). **Known `qiskit_metal` pitfall:** a die-ID label starting with a component name (e.g. `TCZ00` when a component is named `TC`) causes the renderer to silently drop that component's flux-line and route etches. Never start a die ID with a component name.

### Wafer cheese grid and launch clearance profile

For wafer fab exports, [`code/fab/wafer_cheese_grid.py`](../code/fab/wafer_cheese_grid.py) rebuilds ground perforations in the merged positive-metal ground (layer 1/0): the main ground **and every ground island of at least 1 mm²** (e.g. ground enclosed by a CPW loop and tied to the main ground only by air bridges). Smaller metal (signal lines, launches, qubit and coupon pads, all ≤ 0.6 mm² so far) is never perforated. To add cheese to the islands of a finished mask without touching its main ground, use `replace_wafer_cheese(..., islands_only=True)`. Every 5 × 5 µm hole (**40 µm pitch**, 20 µm axis offset so the lattice is symmetric about the die centre, 35 µm metal strip between holes) must fit inside x, y = ±4,452.5 µm (a **547.5 µm unperforated margin** on all four sides of a 10 mm die) and stay ≥ 50 µm from all other etches. The pitch was 20 µm (10 µm offset) until 2026-11-23; masks exported before then keep the 20 µm lattice until re-cheesed. Rationale: a metal strip of width w does not trap vortices when cooled below B ≈ 1.65 Φ₀ / w², so 35 µm strips are safe up to ≈ 2.8 µT residual field (assumes a magnetically shielded fridge; revisit if vortex loss is seen). Standalone test-junction masks keep their own renderer settings with a 500 µm edge keep-out.

Around each populated wafer launch, apply the standard local cheese clearance profile (longitudinal distance measured inward from the nominal PCB anchor to the hole centre; transverse distance from launch centreline to closest hole edge, symmetric on both sides; interpolate between rows and carry the 50 µm row outward to the anchor):

| Inward distance from PCB anchor (µm) | Minimum transverse clearance to a hole edge (µm) |
|---:|---:|
| 50 | 397.5 |
| 100 | 297.5 |
| 150 | 247.5 |
| 200 | 197.5 |
| 250 | 97.5 |
| 300 | 97.5 |

### 10 mm wafer die standard

Every 10 mm wafer die mask follows the same frame:

- **Die boundary & perimeter:** centred 10.000 × 10.000 mm square on layer 1/0, continuous ground perimeter reaching ±5000 µm on both axes.
- **Launches:** nominal PCB anchors from [`pattern design/wafer_pcbport.py`](../pattern%20design/wafer_pcbport.py):
  - Eight-port dies: sides at x = ±4.5 mm, y = ±2.0 mm (P1/P6 at +2.0 mm, P2/P5 at −2.0 mm); top/bottom at x = ±1.1 mm, y = ±4.5 mm.
  - Two-port dies: centred side ports at (±4.5, 0) mm.
  - Spare ports are grounded at their neck.
- **Test-junction coupons:** optional `HookedTestJunction` coupons (F/T styles, Section 6) placed in remaining open ground where they pass `code/layout/check_keepouts.py` (**≥ 100 µm** from both NTHU TiAu and TSRI A12 marker openings and active circuitry); omitted on resonator-only chips without junctions. (The legacy `(±4.20, ±4.00) mm` corner positions conflict with TSRI A12 coarse crosses and are retired for site-portable dies.)
- **Standalone design ID:** placed in the upper-left margin, 50 µm inside both edges (150 µm height, `marker/xmark.py`). Assembled wafers add placement-specific IDs in the upper-right corner ([fabrication_output.md](fabrication_output.md#4-inch-wafer-placement-ids)).
- **Reproducing a wafer die mask:** follow the design-agnostic pipeline in [layout_structure.md §4](layout_structure.md#4-shared-fab-export-pipeline) (`build_design(variant='fab') → export_gds → wafer_cheese_grid.replace_wafer_cheese → check_metal_width.py + check_keepouts.py`). Do not add new Design branches to the legacy `export_wafer_test_junctions.py`.

### Wire-bond launch pad standard

Use the shared `QComponent` in [`pattern design/bondpad.py`](../pattern%20design/bondpad.py): call `add_bond_pad` for nominal PCB anchors (it offsets the taper origin inward along the pad normal), or `BondPad` directly when the taper origin is specified explicitly.

| Parameter | Standard value |
|---|---:|
| Bonding metal width | 400 µm |
| Bonding metal length along the route | 300 µm |
| Side clearance from metal to ground | 200 µm each side |
| Clearance behind the pad | 50 µm |
| Linear taper length | 300 µm |
| Straight lead between taper and `tie` pin | 30 µm |
| Inward shift from nominal PCB anchor | 250 µm |

* The `tie` pin inherits the connected CPW width and gap; the `bond` pin sits on the back edge of the 400 µm pad facing the chip edge across the 50 µm gap.
* On a 10 mm chip with nominal anchors at ±4.5 mm, the etch stops 100 µm from the die edge.
* The 300 µm pad taper is exempt from the generic `< 15°` taper-angle DRC check ([chip_drc.md](chip_drc.md)).
* Apply the local cheese clearance profile above around every populated wafer launch. Preserve existing solver results with their original launch geometry and re-render for any new whole-chip driven validation.

---

## 4. Fab / cheesing settings (`fab=True`)

Implemented in `apply_fab_gds_options()` (`code/layout/gds_tools.py`, design units in mm so `0.005 = 5 µm`, `0.020 = 20 µm`):

| Option | Value | Meaning |
|---|---|---|
| `short_segments_to_not_fillet` | `False` | Fillet every corner (keep every route segment ≥ 2 × fillet). |
| `tolerance` | `1e-6` (1 µm) | Tolerance for filleted path outlines. |
| `fabricate` | `True` | Write only the final mask; drop intermediate keep-out / hole cells. |
| `negative_mask` | `dict(main=[2])` | Layer 1 (all metal) stays **positive**. |
| `no_cheese['buffer']` | `'50um'` | No holes within 50 µm of any etched feature. |
| `no_cheese['view_in_file']['main']` | `{1: False, 2: False}` | Keep-out region not written to file. |
| `cheese['view_in_file']['main']` | `{1: True, 2: False}` | Cheese layer 1 ground only. |
| `cheese_0_x`, `cheese_0_y` | `0.005` (5 × 5 µm) | Square holes. |
| `delta_x`, `delta_y` | `0.020` (20 µm) | Legacy standalone renderer pitch; **wafer fab masks rebuild a 40 µm lattice** via `code/fab/wafer_cheese_grid.py` (Section 3). |
| `cheese['edge_nocheese']` | `'0.5mm'` | Base renderer margin (wafer exports enforce 547.5 µm via `wafer_cheese_grid.py`). |

Fab GDS layers: `1/100` cheesed ground (`TOP_main_1_Cheese_100`), `1/10` and `1/11` pads/islands/CPW centre conductors. Without cheesing, ground is on `1/0`.

---

## 5. Units and cell cleanup

`qiskit_metal` writes mm; `_prune_gds` re-reads with `units='convert'` and writes **µm** (`unit=1e-6`, `precision=1e-9`). (Using `gdspy`'s default `units='skip'` keeps mm numbers labelled as µm, making a 10 mm chip appear 10 µm wide in KLayout.) Unreferenced cells and the dangling `Cheese_diff` reference are pruned.

---

## 6. Test-junction coupons

<a id="6-corner-test-junction-coupons"></a>

Chips with Josephson junctions can include probe coupons for fabrication process monitoring wherever open ground remains (resonator-only chips without junctions omit them):
- **Component:** `HookedTestJunction` ([`pattern design/junction_test_coupon.py`](../pattern%20design/junction_test_coupon.py)), declared directly in each Design's `placement.json` ([layout_structure.md](layout_structure.md)).
- **Geometry:** two floating 100 × 300 µm pads with 10 µm landing stems, 2 × 4 µm Dolan recesses in each stem tip when `fab_contacts=True`, inside an etched 230 × 820 µm pocket.
- **Styles:**
  - **F** (fluxonium): 30 µm opening with recesses facing each other (room for a 60 µm array envelope).
  - **T** (transmon): 5 µm opening with the upper recess offset by +2 µm.
- **Placement & clearance:** place coupons in open ground regions so every pocket passes `code/layout/check_keepouts.py` (**≥ 100 µm** from both NTHU TiAu and TSRI A12 marker openings and from active circuitry/routes). Avoid the legacy `(±4.20, ±4.00) mm` corner slots, which collide with TSRI A12 coarse crosses.
- **Junction:** never drawn on the base-metal mask; 50 µm cheese keep-out applies around the pocket.
