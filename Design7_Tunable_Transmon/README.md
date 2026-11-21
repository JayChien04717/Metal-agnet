# Design7_Tunable_Transmon: 8-port half-resonator / half-tunable-transmon chip

> **Status:** Stage 1 layout draft (`2026-11-20_r1`, unsimulated). Waiting for visual layout review before running EM simulations (Q3D / HFSS) or generating final fabrication masks.

**Figures:** [fig/](fig/README.md) holds the whole chip ([`chip.png`](fig/chip.png)) and part close-ups ([`bare_resonator.png`](fig/bare_resonator.png), [`q1.png`](fig/q1.png), [`q1_pocket.png`](fig/q1_pocket.png), [`q1_squid.png`](fig/q1_squid.png)), rendered from the companion simulation GDS by `code/general/plot_design_figs.py`.

## 1. Design purpose and topology

A 10 × 10 mm 8-port chip (`wafer_pcb_ports`, side ports at `y = ±2.0 mm`) split into two independent feedlines:

- **Top half (`P1 ↔ P6`, `y = +2.0 mm`) — Resonator-only:** 4 bare λ/4 hanger resonators (`R1`–`R4`, analytic targets `6.40`–`7.00 GHz`, `10 / 6 µm` CPW, `300 µm` coupling length, `8 µm` ground strip to feedline).
- **Bottom half (`P2 ↔ P5`, `y = -2.0 mm`) — Tunable transmons (no charge lines, no flux lines):** 4 custom `TunableTransmonPocket` qubits (`Q1`–`Q4`) multiplexed on 4 λ/4 readout hanger resonators (`RR1`–`RR4`, analytic targets `6.50`–`7.10 GHz`). Flux tuning is external (coil/magnet); readout and drive go through the feedline.
- **Spare ports (`P3`, `P4`, `P7`, `P8`):** standard wire-bond launch pads grounded at their neck.

## 2. Custom component (`TunableTransmonPocket`)

Defined locally in [`layout/components/tunable_transmon_pocket.py`](layout/components/tunable_transmon_pocket.py) (candidate for promotion to `pattern design/transmon/` after user verification):

- Floating double-pad transmon pocket (`450 × 120 µm` pads, `50 µm` pad gap, `15 µm` pad corner radius, `650 × 550 µm` rounded pocket) with optional top/bottom readout coupling pad (`connection_pads`).
- Two vertical SQUID landing legs (`16 µm` wide, `30 µm` center-to-center `squid_spacing`, `5 µm` `jj_cut`) bridging `pad_top` and `pad_bot`.
- `TunableTransmonPocketFab` (`fab_contacts='True'`) cuts `2 × 4 µm` Dolan contact hooks (lower hooks offset by `+2 µm` in x) while maintaining ≥ 5 µm metal on both sides of every hook (passing the 4 µm minimum metal-width DRC C12a).
- Implements `contact_hooks()` directly compatible with `code/general/ebeam_junction_placement.py` (`add_ebeam_junctions()`) and `pattern design/ebeam_junction/transmon_junction.py` (`TransmonJunction`). All 4 qubits keep `orientation='0'` so every Dolan bridge shares the vertical finger-left evaporation direction.

## 3. Artifacts and checks (`2026-11-20_r1`)

- Simulation companion GDS: [`gds/Design7_2026-11-20_r1_sim.gds`](gds/Design7_2026-11-20_r1_sim.gds) (lumped junction annotations on layer `50/0`, unnotched).
- E-beam junction placement preview GDS: [`gds/Design7_2026-11-20_r1_ebj_preview.gds`](gds/Design7_2026-11-20_r1_ebj_preview.gds) (`TunableTransmonPocketFab` hooks + placeholder `TransmonJunction` on layers `11/0` and `12/0`) · [Q1 overlay check](notes/fig/q1_transmon_junction_overlay.png).
- **Layout DRC checks passed:**
  - All 8 filleted resonator lengths match analytic targets to `< 0.001 µm`.
  - Every route bend uses `fillet = 70 µm` with all straight segments `≥ 180 µm > 2 × fillet`.
  - All 4 qubit `readout` pin normals point outward cleanly.
  - Windowed 4 µm minimum metal-width check (`code/general/check_metal_width.py`): `0` violations on both simulation and fab-hooks GDS.
  - **No EM simulations (Q3D / HFSS) have been run yet.**

## 4. Parameter files and reproduction

- [design.json](design.json) · [geometry.json](parameters/geometry.json) · [quantum.json](parameters/quantum.json) · [solver.json](parameters/solver.json) · [provenance.json](parameters/provenance.json)

Reproduce from the repository root (`qmetal` environment):

```powershell
python -X utf8 -B Design7_Tunable_Transmon/layout/build.py --ebeam-preview
python -X utf8 -B code/general/plot_design_figs.py Design7_Tunable_Transmon
```
