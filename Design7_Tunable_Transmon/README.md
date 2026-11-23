# Design7_Tunable_Transmon: 8-port half-resonator / half-tunable-transmon chip

> **Status:** Stage 2 single-qubit Q3D and single-channel HFSS driven S21 tuned (`2026-11-22_r2`). Full-chip multi-channel EM validation and final fabrication GDS confirmation remain pending.

**Figures:** [fig/](fig/README.md) holds the whole chip ([`chip.png`](fig/chip.png)) and part close-ups ([`bare_resonator.png`](fig/bare_resonator.png), [`q1.png`](fig/q1.png), [`q1_pocket.png`](fig/q1_pocket.png), [`q1_squid.png`](fig/q1_squid.png)), rendered from the companion simulation GDS by `code/general/plot_design_figs.py`.

![Design7_Tunable_Transmon whole chip (2026-11-22_r2)](fig/chip.png)

| `TunableTransmonPocket` Q1 close-up | Bare λ/4 hanger resonator R1 |
|---|---|
| ![Q1 pocket](fig/q1_pocket.png) | ![Bare resonator R1](fig/bare_resonator.png) |

## 1. Design purpose and topology

A 10 × 10 mm 8-port chip (`wafer_pcb_ports`, side ports at `y = ±2.0 mm`) split into two independent straight feedlines with alternating (`+1, -1, +1, -1`) channel placement so bare resonators and qubits are interleaved horizontally in the central ground region:

- **Top half (`P1 ↔ P6`, `y = +2.0 mm`) — Resonator-only:** 4 bare λ/4 hanger resonators (`R1`–`R4`, targets `6.30, 6.70, 6.50, 6.90 GHz` — same frequency set as the bottom half, `10 / 6 µm` CPW, `300 µm` coupling length, `4 µm` ground strip to feedline).
- **Bottom half (`P2 ↔ P5`, `y = -2.0 mm`) — Tunable transmons (no charge lines, no flux lines):** 4 custom `TunableTransmonPocket` qubits (`Q1`–`Q4`, symmetric SQUIDs, `orientation='0'`) multiplexed on 4 λ/4 readout hanger resonators (`RR1`–`RR4`, targets `6.30, 6.70, 6.50, 6.90 GHz`). Flux tuning is external (coil/magnet); drive and readout go through the feedline.
- **Spare ports (`P3`, `P4`, `P7`, `P8`):** standard wire-bond launch pads grounded at their neck.

## 2. Parameter summary and provenance (`2026-11-22_r2`)

### Qubit & dispersive readout (Q3D + exact transmon diagonalization)

Source: single-qubit Q3D run [`20261122_201149_q3d_ec200_v2_34009f36`](results/qubit/20261122_201149_q3d_ec200_v2_34009f36/summaries/summary.json) (8 passes, 0.5% target, Medium order, `Q1` + open-terminated readout pad). `EJ` and `Lj` are quantum-model design inputs (not HFSS-extracted quantities); e-beam mask widths are estimates from `code/qubit/jj_resistance_model.py` (80 nm finger, $R_A = 625\,\Omega\cdot\mu\text{m}^2$).

| Parameter | Target | Adopted value | Unit | Role / status |
|---|---:|---:|---|---|
| Charging energy $E_C / h$ | 200 | **199.68** | MHz | Result (Q3D, $C_\Sigma = 97.00\text{ fF}$) |
| Qubit frequency $f_{01}$ ($\Phi = 0$) | ~4.00 | **4.000** | GHz | Derived from Q3D $E_C$ + input $E_J$ |
| Anharmonicity $\alpha / 2\pi$ | — | **−227.4** | MHz | Derived (exact charge-basis diag.) |
| Symmetric SQUID $E_{J,\mathrm{sum}} / h$ | — | **11.098** | GHz | Input ($E_{J1} = E_{J2} = 5.549\text{ GHz}$) |
| Branch inductance $L_{J1} = L_{J2}$ | — | **29.46** | nH | Input (HFSS lumped branch value) |
| Predicted branch $R_n$ ($T \to 0$) | — | **25.31** | kΩ | Derived (`jj_resistance_model.py`) |
| `TransmonJunction` widths (finger / cross) | — | **80 / 155** | nm | Derived (symmetric: `jj1 = jj2 = 155 nm`) |
| Readout voltage-division ratio $\beta$ | — | **0.07069** | — | Result (Q3D, $C_{\mathrm{ro,top}} = 15.69\text{ fF}$) |
| Readout pad self-capacitance $C_{\mathrm{ro,self}}$ | — | **35.98** | fF | Result (Q3D → `0.220 mm` length shortening) |
| Coupling $g / 2\pi$ (`6.3–6.9 GHz`, $\lambda/4$) | — | **88.5–96.9** | MHz | Derived from Q3D $\beta$ and $\lambda/4$ $V_{\mathrm{zpf}}$ |
| Dispersive shift $|2\chi| / 2\pi$ (`6.3–6.9 GHz`) | 0.50 | **0.471–0.613** | MHz | Derived (`0.556` at 6.5 GHz, `0.510` at 6.7 GHz) |

| Channel | $f_{r,\mathrm{target}}$ (GHz) | $L_{\mathrm{res}}$ (mm) | $g / 2\pi$ (MHz) | $2\chi / 2\pi$ (MHz) |
|---|---:|---:|---:|---:|
| `Q1 / RR1` | 6.30 | 4.4379 | 88.49 | −0.613 |
| `Q3 / RR3` | 6.50 | 4.2946 | 91.30 | −0.556 |
| `Q2 / RR2` | 6.70 | 4.1599 | 94.11 | −0.510 |
| `Q4 / RR4` | 6.90 | 4.0329 | 96.92 | −0.471 |

### Hanger resonators (HFSS driven S21)

Source: single-channel HFSS driven S21 run [`20261122_205320_hanger_c300_s4_verify_35cbb6d1`](results/resonator/20261122_205320_hanger_c300_s4_verify_35cbb6d1/summaries/summary.json) (4 adaptive passes, $\Delta S = 0.01$ setup, explicit `10 µm` length-based mesh on `trace_feedline` and `trace_R_test`).

| Parameter | Target | Simulated value | Unit | Role / status |
|---|---:|---:|---|---|
| Bare calibration $f_r$ (`L = 4.5146 mm`) | 6.500 | **6.49836** | GHz | Result (HFSS driven S21, −1.64 MHz) |
| Calibrated phase velocity $v_{\mathrm{eff}}$ | — | **1.1738 × 10⁸** | m/s | Derived from HFSS $f_r$ |
| External quality factor $Q_c$ (`300 µm / 4 µm`) | 13,000 | **12,346** | — | Result (HFSS complex notch fit, $Q_l = 12,346$) |
| Linewidth $\kappa / 2\pi$ at 6.50 GHz | 0.50 | **0.526** | MHz | Result (HFSS driven S21) |

| Bare channel | $x$ (mm) | Side | $f_{\mathrm{target}}$ (GHz) | Filleted $L$ (mm) |
|---|---:|---:|---:|---:|
| `R1` | −2.55 | +1 (up) | 6.30 | 4.6579 |
| `R2` | −0.85 | −1 (down) | 6.70 | 4.3799 |
| `R3` | +0.85 | +1 (up) | 6.50 | 4.5146 |
| `R4` | +2.55 | −1 (down) | 6.90 | 4.2529 |

### Geometry knobs (`parameters/geometry.json`)

- **TunableTransmonPocket:** `pad_width = 595 µm`, `pad_height = 165 µm`, `pad_gap = 38 µm`, `pad_corner_radius = 15 µm`, `pocket_width = 795 µm`, `pocket_height = 620 µm`, `pocket_corner_radius = 25 µm`.
- **SQUID landing legs:** `squid_width = 46 µm`, `squid_loop_width = 14 µm` (`16 µm` legs, `30 µm` `squid_spacing`, `5 µm` `jj_cut`); `TunableTransmonPocketFab` adds `2 × 4 µm` Dolan hooks (`+2 µm` lower x-offset) compatible with `add_ebeam_junctions()` and `TransmonJunction`.
- **Readout coupling pad:** `165 × 10 µm` T-bar, `14 µm` gap to transmon pad, `50 µm` CPW extension.
- **Hanger resonator:** `10 / 6 µm` CPW, `300 µm` coupling length, `4 µm` coupling ground strip, `300 µm` meander pitch (`278 µm` trench gap), `70 µm` fillet.

## 3. Artifacts, DRC, and solver overrides

- Simulation companion GDS: [`gds/Design7_2026-11-22_r2_sim.gds`](gds/Design7_2026-11-22_r2_sim.gds) (lumped junction sheets on layer `50/0`).
- E-beam junction placement preview GDS: [`gds/Design7_2026-11-22_r2_ebj_preview.gds`](gds/Design7_2026-11-22_r2_ebj_preview.gds) (`TunableTransmonPocketFab` hooks + `TransmonJunction` `80 / 155 nm` on layers `11/0` and `12/0`).
- **DRC checks passed:** all 8 resonator filleted lengths match target to `< 0.001 µm`; all straight segments `≥ 150 µm > 2 × 70 µm` fillet; all 4 qubit `readout` pins point outward; windowed 4 µm minimum metal-width check (`code/general/check_metal_width.py`): `0` violations on both simulation and fab-hooks GDS.
- **Design-specific solver overrides:**
  - AEDT project filename is shortened to `simulation/aedt/D7.aedt` so `.aedtresults` paths remain below the Windows 260-character `MAX_PATH` limit.
  - Auto wirebonds (`wb_size = 2`) are disabled on the single-channel hanger coupling fixture (`hfss_hanger.py`) because `wb_size = 2` bridges across the `4 µm` coupling ground strip and shorts the feedline to the resonator.
- **Pending validation:** full-chip loaded readout resonator HFSS check (to confirm exact `0.220 mm` capacitive-load shortening on `RR1`–`RR4`), Q3D pass-count convergence check beyond 8 passes, and Stage 3 fab GDS checklist confirmation.

## 4. Parameter files and reproduction

- [design.json](design.json) · [geometry.json](parameters/geometry.json) · [quantum.json](parameters/quantum.json) · [solver.json](parameters/solver.json) · [provenance.json](parameters/provenance.json)

Reproduce from the repository root (`qmetal` environment):

```powershell
# Layout, GDS, and figures
python -X utf8 -B Design7_Tunable_Transmon/layout/build.py --ebeam-preview
python -X utf8 -B code/general/plot_design_figs.py Design7_Tunable_Transmon

# Single-qubit Q3D (AEDT running)
python -X utf8 -B Design7_Tunable_Transmon/simulation/qubit/q3d_qubit.py --tag q3d_ec200_v2

# Single-channel HFSS driven S21 (AEDT running)
python -X utf8 -B Design7_Tunable_Transmon/simulation/resonator/hfss_hanger.py --tag hanger_c300_s4_verify --f-target 6.5 --coup-len 300 --coup-strip 4 --v-ph 1.1738e8
```
