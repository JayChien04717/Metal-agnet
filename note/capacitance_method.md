# Capacitance-matrix method for one FTF set

Code: [cap_matrix.py](cap_matrix.py). Simulation: [../simulation/q3d_ftf.py](../simulation/q3d_ftf.py).

## 1. Q3D run

The simulation renders FX1, TC and FX2 (`build_design(for_process=True)`, so each
fluxonium's JJ and array are replaced by one lumped inductor). Every pin is given
an open termination, so the readout pads and charge lines become separate
conductors. Q3D returns the Maxwell capacitance matrix **C** in fF.

Net names are `<polygon>_<component>`, for example `pad_top_FX1`, `tmon_TC`,
`readout_line_FX1` and `ground_main_plane`.

## 2. Nodes

| index | node | meaning |
|---|---|---|
| 0, 1 | `pad_bot_FX1`, `pad_top_FX1` | FX1 pads → sum / difference |
| 2 | `tmon_TC` | coupler island (grounded transmon) |
| 3, 4 | `pad_top_FX2`, `pad_bot_FX2` | FX2 pads → sum / difference |
| 5 | `readout_line_FX1` | resonator 1 |
| 6 | `readout_line_FX2` | resonator 2 |

Any conductor left out of this list is grounded: its rows and columns are dropped.
That covers ground, the flux lines, and the charge lines, which are connected to
50 Ω at the qubit frequency. The thesis uses an 8 × 8 matrix. This design has no
coupler readout, so it uses 7 × 7; the method is the same.

## 3. Sum / difference coordinates (thesis A.14)

M has the block [[1, 1], [1, -1]] on each FX pad pair and 1 elsewhere on the diagonal.

    C~ = M^-1ᵀ · C · M^-1

## 4. Read-off (thesis A.15, A.18), 0-based indices

    Ec_FX1 = e²/2 · C~^-1[1,1]      J_1c = 4e² · C~^-1[1,2]
    Ec_TC  = e²/2 · C~^-1[2,2]      J_2c = 4e² · C~^-1[2,4]
    Ec_FX2 = e²/2 · C~^-1[4,4]      J_12 = 4e² · C~^-1[1,4]

    ħ g_FX1-R1 = e · sqrt(π ħ / 2 Z0) · C~^-1[1,5]      (Z0 = 50 Ω)
    ħ g_FX2-R2 = e · sqrt(π ħ / 2 Z0) · C~^-1[4,6]

For a fluxonium, g excludes the qubit matrix element (A.17). For a transmon,
multiply by (32 Ec / EJ)^(1/4), which is about 0.9 at EJ/Ec ≈ 50 (A.21).

## 5. Targets

| quantity | target | total capacitance C_Σ = e² / (2 h Ec) |
|---|---|---|
| Ec_FX | 1.3 GHz | 14.9 fF |
| Ec_TC | 0.2 GHz | 96.9 fF |

## 6. Results

Every run is in `../simulation/results/` (see `simulation/README.md`): the runs behind the
quoted numbers in `final/`, earlier ones in `history/<stage>/`, as `<tag>_cmat_fF.csv` (Maxwell matrix)
plus `<tag>_summary.json`, and each sweep has a `sweep_<step>.csv` table.

Q3D settings: Medium solution order, no auto increase, 0.5 % error, at most 6
adaptive passes. Each run takes about 0.5–1 min.

### Current design (v4, `v4_final`): values to quote

Changes in v4:
- **FX pocket:** shrunk to 260 × 360 µm and made asymmetric (`pocket_top` 210 µm,
  `pocket_bottom` 150 µm), so the bottom edge is 40 µm from the bottom pad.
- **FX charge line:** back behind a ground strip (`cl_sep` 15 µm, ~9 µm of ground
  after its end gap). The closer pocket edge keeps the coupling strong.
- **FX readout pad:** width halved to 100 µm. The pad-to-pad gap is set per qubit
  (`FX_PER_QUBIT`: 19.75 µm for FX1, 21.2 µm for FX2) so that g = 150 MHz for both
  resonator frequencies.
- **FX pads:** 125 × 50 µm.
- **Coupler:** stem 158 µm, readout pad 273 µm.
- **Resonator lengths:** recomputed with the measured readout-pad capacitances
  (`v4_pads_only`: 24.5 / 24.1 / 48.8 fF for FX1 / FX2 / TC) → 9.302 / 9.015 / 9.606 mm.

Sweeps: `results/history/6_v4_pocket_readout/sweep_fx1.csv` and runs `v4_base`, `v4_tuned`,
`v4_final` (final run in `results/final/`).

| quantity | value | target |
|---|---|---|
| Ec_FX1 / Ec_FX2 | 1.309 / 1.313 GHz | 1.3 GHz |
| Ec_TC | 0.202 GHz | 0.2 GHz |
| J_1c / J_2c / J_12 | 138 / 139 / 12 MHz | – |
| g/2π FX1–R1 / FX2–R2 (per n, A.17) | 150.1 / 148.9 MHz (β = 0.202 / 0.194) | 150 MHz |
| g/2π TC–R_TC (transmon convention, EJ/Ec = 50) | 148.6 MHz | 150 MHz |
| κ/2π | 2 / 3 / 2 MHz (by design, tee fit) | 2 / 3 / 2 MHz |
| FX charge line | 0.45 fF, β 0.0129, T1 limit 2.0 ms, −86 dBm for 10 MHz Rabi | – |
| TC charge line | 0.047 fF, T1 limit 1.1 ms | – |

Residual deviations of ≤ 1 % are at the level of the 6-pass Q3D mesh noise.

### v3 (`v3_final` / `v3b_final` with the stronger FX charge line)

v3 adds the coupler readout and three coupled-line tees, with resonators at
6.3 / 6.0 / 6.5 GHz for FX1 / TC / FX2. All three resonators are rendered in Q3D.
The coupler was re-tuned (sweeps `results/history/4_v3_resonators_TC_readout/sweep_tc1..tc3.csv`): stem 250 → 170 µm, readout
pad 272 × 20 µm placed 5 µm from the bar, charge-line tip 250 µm below the bar etch.

| quantity | value | target |
|---|---|---|
| Ec_FX1 / Ec_FX2 | 1.315 / 1.312 GHz | 1.3 GHz |
| Ec_TC | 0.1986 GHz | 0.2 GHz |
| J_1c / J_2c / J_12 | 124 / 123 / 11 MHz | – |
| g/2π FX1–R1 / FX2–R2 (per n, A.17 convention) | 167 / 173 MHz | – |
| g/2π TC–R_TC, transmon convention (A.21, EJ/Ec = 50) | 151 MHz | 150 MHz |
| κ/2π (tee Cc from Q3D fit, section 8) | 2 / 3 / 2 MHz (by design) | 2 / 3 / 2 MHz |

Charge lines (`cap_matrix.charge_line`; 50 Ω line, T = 0; drive power for a 10 MHz Rabi rate):

| line | C to qubit | β | T1 limit | P at chip for 10 MHz Rabi |
|---|---|---|---|---|
| FX1 (f01 0.5 GHz, \|n01\| 0.1 assumed) | 0.50 fF | 0.0141 | 1.6 ms | −87 dBm |
| TC (5 GHz, EJ/Ec 50 assumed) | 0.049 fF | 0.0005 | 1.06 ms | −79 dBm |

With the tip 15 µm or 80 µm from the bar, the TC charge line gave 0.75 / 0.38 fF
(T1 limit 4 / 21 µs), which is too strong.

FX charge line: the first version (tip 15 µm outside the pocket, 0.18 fF, T1 limit
14 ms) needed −77 dBm at the chip. It was strengthened about 3× by moving the tip
14 µm *into* the pocket (`charge_line_options.cl_sep = '-14um'`, run `v3b_final`).
That saves ~10 dB of drive power and keeps the T1 limit at ~1.6 ms. Sweep
(`simulation/tuning/sweep_FX_charge_line.py`, `results/history/5_v3_FX_charge_line/sweep_fx_charge.csv`):

| cl_sep | C to FX | β | T1 limit | P for 10 MHz Rabi |
|---|---|---|---|---|
| +15 µm (v3) | 0.18 fF | 0.0047 | 14.5 ms | −77 dBm |
| +5 µm | 0.29 fF | 0.0080 | 5.1 ms | −82 dBm |
| −10 µm | 0.45 fF | 0.0128 | 2.0 ms | −86 dBm |
| **−14 µm (chosen)** | **0.50 fF** | **0.0141** | **1.6 ms** | **−87 dBm** |
| −30 µm | 0.72 fF | 0.0210 | 0.74 ms | −90 dBm |

The T1 limit and drive power assume f01 = 0.5 GHz and |n01| = 0.1. They scale as
T1 ∝ 1 / (|n01|² f01) and P ∝ 1 / |n01|², so redo them with the real EJ, EL.
`v3b_final` leaves the other numbers unchanged (Ec_FX 1.313 / 1.310 GHz,
Ec_TC 0.1985 GHz, g_TC 151 MHz).

### v2 (`v2_with_res`)

Changes in v2: JJ at the pocket centre, 50 µm hairpin array fully inside the
pocket (no array slot, so no ground near the array), and a U-loop flux line from the
pocket's left wall with its tip 20 µm from the array's end. Resonators are included
in Q3D, and g uses the design f_r.

| quantity | value | target |
|---|---|---|
| Ec_FX1 / Ec_FX2 | 1.308 / 1.311 GHz | 1.3 GHz |
| Ec_TC | 0.2015 GHz | 0.2 GHz |
| J_1c / J_2c | 122 / 124 MHz | – |
| J_12 | 12.5 MHz | – |
| g/2π FX1–R1 / FX2–R2 | 187 / 195 MHz (β = 0.227 / 0.225) | – |

The flux line's U is shorted to ground at the wall, so it is the nearest
ground-connected metal to the array, 20 µm from the array's end
(`flux_line_options.fl_sep`). The pocket walls are ≥ 90 µm away.

### v1, pads only (`final_default`, array in a slot); see section 7 for the values with resonators

| quantity | value | target |
|---|---|---|
| Ec_FX1 / Ec_FX2 | 1.319 / 1.320 GHz | 1.3 GHz |
| Ec_TC | 0.1997 GHz | 0.2 GHz |
| J_1c / J_2c | 157 / 157 MHz | – |
| J_12 | 17 MHz | – |

Geometry: FX pads 130 × 50 µm, pad gap 120 µm, FX pocket 300 × 480 µm, FX
centres at x = ±0.34 mm. Coupler bar 15 µm wide (arm length 370 µm), stem 250 µm.

### Tuning log

| run | change | Ec_FX [GHz] | Ec_TC [GHz] | C_Σ FX / TC [fF] |
|---|---|---|---|---|
| baseline | pads 240×90, FX_X 0.55 | 0.648 | 0.112 | 29.9 / 173 |
| s1_pad160x60 | pads 160×60 | 1.044 | 0.111 | 18.6 / 174 |
| s1_pad120x50 | pads 120×50 | 1.350 | 0.110 | 14.3 / 176 |
| s1_pad120x50_ro | + smaller readout pad | 1.356 | 0.110 | 14.3 / 176 |
| s1_tc_short | stem 250, bar 15 µm | 0.649 | 0.134 | 29.9 / 145 |
| s2_x035 | pads 130×50, pocket 300, FX_X 0.35 | 1.312 | 0.195 | 14.8 / 99.4 |
| s2_x030 | same, FX_X 0.30 | 1.335 | 0.225 | 14.5 / 86.1 |
| s3_x034 | same, FX_X 0.34 | 1.312 / 1.320 | 0.200 | 14.8 / 96.9 |

What drives each Ec:
- FX Ec is set almost entirely by pad size.
- Coupler Ec is set by the bar length, i.e. how far apart the two FXs are.
- The readout pad size barely changes Ec.

## 7. Readout coupling g, with the resonators in Q3D (`final_with_res`, v1 geometry; v2 numbers are in section 6)

Command: `Q3DRunner().run(tag, with_resonators=True)`. This renders `res1` and `res2`
as well. Their tee ends are open; at the qubit end each resonator merges with the
readout pad into one net (`readout_line_FX1/2`, C_line ≈ 1.19 / 1.12 pF).

Method (`cap_matrix.readout_g`):

    β      = C~^-1[q, r] / C~^-1[r, r]                  (= C_c / C_Σ,q, voltage division)
    V_zpf  = ω_r · sqrt(ħ Z0 / π)                     (λ/2 resonator, open end)
    ħ g    = 2e · β · V_zpf                           (per charge number n, as in A.17)

For a lumped LC this is identical to A.18, ħg = e·sqrt(2ħ/Z)·C~^-1[q, r]. A.18 itself
is written for a λ/4 resonator (Z = 4Z0/π) and needs node r to carry the
resonator's *effective* capacitance. Plugging the full DC capacitance from Q3D into
A.18 gives 84 / 89 MHz, about 2.2× too low for our λ/2 resonators.

f_r comes from the design: FXTFX24's measured f·L ≈ 53 GHz·mm with RES_LEN = 7.6 / 7.2 mm
gives 6.97 / 7.36 GHz. The estimate 1/(2·Z0·C_line) from Q3D gives 9.3–9.9 GHz, which
is too high because the tee section isn't rendered and the 6-pass mesh is coarse.

| | β | f_r [GHz] | V_zpf [µV] | g/2π [MHz] |
|---|---|---|---|---|
| FX1 – R1 | 0.224 | 6.97 | 1.70 | 184 |
| FX2 – R2 | 0.224 | 7.36 | 1.80 | 194 |

Multiply g by the fluxonium charge matrix element of the transition you read out
to get the actual coupling. To make g smaller, move the readout pad away from the
top pad (`readout_line_options.pad_sep`) or shrink it.

With the resonator attached, the readout node is no longer a small floating pad,
so the other numbers become (`final_with_res`, also in `results/history/2_v1_readout_g/final_with_res_g.json`):

| Ec_FX1 / Ec_FX2 | Ec_TC | J_1c / J_2c | J_12 |
|---|---|---|---|
| 1.304 / 1.304 GHz | 0.201 GHz | 124 / 125 MHz | 15 MHz |

These are the values to quote. The pad-only run (`final_default`) overestimates J
through the floating readout pads.

## 8. Resonator linewidth κ (coupling to the feedline)

Lab formula, from FXTFX24-W7-0x-V.3.1 cell 42:

```python
# C_1   : coupling capacitance resonator <-> feedline from Q3D (F)
# Z_0   : feedline impedance, Z_r : resonator impedance (both 45 Ohm in this design)
omega_1 = 2*pi*f_1
Q_c1    = pi / (Z_0 * Z_r * (C_1 * omega_1)**2)      # coupling quality factor
Kappa_1 = (f_1 / Q_c1) * 1e-6                       # linewidth kappa/2pi in MHz
```

Inverted, to design for a target κ: `Q_c = f / κ`, `C_c = sqrt(π / (Q_c Z_0 Z_r)) / ω`.
This is `coupling_cap()` in `pattern design/ftf_layout.py`.

How it is used here:

1. **Cc vs tee coupling length.** `simulation/q3d_tee.py` runs Q3D on one
   CoupledLineTee (15 / 5.909 µm CPW, 3 µm ground strip, 150 µm drop). Cc is the
   mutual capacitance between `prime_cpw` (feedline) and `second_cpw` (resonator).
   Results are in `results/final/tee_sweep.csv`:

   | coupling length | 0.4 mm | 0.8 mm | 1.2 mm |
   |---|---|---|---|
   | Cc | 8.49 fF | 18.16 fF | 27.67 fF |

   Linear fit: **Cc = 23.985 fF/mm · L − 1.082 fF** (`CC_SLOPE`, `CC_OFFSET` in ftf_layout.py).

2. **Coupling length for each target κ** (Z_0 = Z_r = 45 Ω):

   | resonator | f | κ/2π target | Q_c | Cc | coupling length |
   |---|---|---|---|---|---|
   | FX1 | 6.3 GHz | 2 MHz | 3150 | 17.73 fF | 0.784 mm |
   | FX2 | 6.5 GHz | 2 MHz | 3250 | 16.92 fF | 0.750 mm |
   | TC  | 6.0 GHz | 3 MHz | 2000 | 23.36 fF | 1.019 mm |

3. **Resonator length.** λ/2, open at both ends, from the analytic CPW
   (Z0 = 45.1 Ω, C = 184.5 pF/m, v = 1.2016e8 m/s):
   f = v / 2L · (1 − C_load / C_line), with C_load = readout pad + Cc. The pad
   capacitances are measured in Q3D (`v4_pads_only`: 24.5 / 24.1 / 48.8 fF). The
   length is measured along the real filleted centre line, from the tee's open end
   to the qubit pad, and the routes are iterated in `build_design` until it hits the
   target: 9.302 / 9.606 / 9.015 mm (FX1 / TC / FX2).

### 8b. κ from the capacitance matrix of the full readout circuit

`simulation/q3d_kappa.py` renders the feedline, all three tees and all three
resonators together in Q3D; the qubits are left out and the resonators' qubit ends
are open. C_1 = −C[resonator net, feedline net] is read from the Maxwell matrix and
put into the lab formula above. This is how FXTFX24 did it (cell 41: TL + cavities).
It also captures the coupling of the rest of each resonator to the feedline, which
the isolated-tee fit misses.

- **Mesh:** κ ∝ C_1², and this structure is large. With 6 passes C_1 scattered by
  ±0.5 fF (±6 % in κ), so κ runs use 12 passes at 0.2 % (`KAPPA_SETUP`, ~2.5 min).
  Even then the run-to-run noise is about ±0.3 fF, or ±3 % in κ.
- **Correction:** the full-circuit C_1 minus the isolated-tee fit is stored per
  resonator as `C_extra` in `RESONATORS`. The tee length is then
  (Cc_target − C_extra − CC_OFFSET) / CC_SLOPE.

| run | tee FX1 / TC / FX2 [mm] | C_1 [fF] | κ/2π [MHz] |
|---|---|---|---|
| kappa_v4 (6 passes) | 0.784 / 1.019 / 0.750 | 18.64 / 25.30 / 17.51 | 2.21 / 3.52 / 2.14 |
| kappa_v4c_conv (12 passes) | 0.749 / 0.979 / 0.710 | 17.11 / 22.98 / 16.12 | 1.86 / 2.90 / 1.82 |
| **kappa_v4d_conv (12 passes, final)** | **0.775 / 0.995 / 0.743** | **18.10 / 23.56 / 17.23** | **2.08 / 3.05 / 2.08** |

The final κ is within 2–4 % of the targets (2 / 3 / 2 MHz), about the size of the
mesh noise. The tee changes don't touch the qubits, so the v4 Ec / g values stand.

Not yet checked: f_r has not been verified with an HFSS eigenmode run. The
1 − C_load/C_line loading correction and the lab κ formula are first-order estimates.

## 9. Resonator frequency from HFSS eigenmode

Script: `simulation/hfss_eigen_res.py` (qiskit-metal `EigenmodeSim`, HFSS 2021 R2).

- **What is rendered, per resonator:** the resonator route, its coupled-line tee
  (feedline section open at both ends, so its own modes are > 9 GHz), and its qubit.
  The qubit is the simulation version with lumped inductors, so the readout pad
  loads the resonator as in the device.
- **Qubits kept far below the resonators:**
  - Fluxonium: lumped 200 nH, ~3 GHz.
  - Coupler: `Lj` = 50 nH per junction (new `MyTmonCoupler` option `Lj`), ~2.3 GHz.

  The resonator mode is then only dispersively pulled upward, by ~g²/Δ ≈ 5–10 MHz.
  The real readout frequency also depends on where the qubits are operated.
- **Setup:** `min_freq_ghz` 4, 2 modes, `max_delta_f` 0.05 %, at most 17 passes,
  30 % refinement. About 8–11 minutes per resonator. 25 passes pushed `hf3d` above
  14 GB of RAM, so stay at ≤ 17 on this PC. The last passes still move ~3–5 MHz each,
  so treat results as ±10 MHz.

| resonator | length (analytic) | f target | f HFSS | error |
|---|---|---|---|---|
| FX1 | 9.302 mm | 6.300 GHz | 6.2892 GHz | −10.8 MHz (−0.17 %) |
| TC  | 9.606 mm | 6.000 GHz | 6.0259 GHz | +25.9 MHz (+0.43 %) |
| FX2 | 9.015 mm | 6.500 GHz | 6.4892 GHz | −10.8 MHz (−0.17 %) |

Second mode in each run (9.0–9.3 GHz): from the open tee / feedline section, not a resonator harmonic.

**Correction:** f ∝ 1/L, so each length is multiplied by f_HFSS / f_target
(`hfss_scale` in `RESONATORS`). New lengths: **9.286 / 9.647 / 9.000 mm**
(FX1 −16 µm, TC +41 µm, FX2 −15 µm). The corrected lengths were not re-simulated.
The correction is a ≤ 0.4 % rescale, so the expected residual is within the ±10 MHz
mesh uncertainty. Fabrication spread (kinetic inductance, ε_r, etch) is typically
larger than that, so leave some margin in the readout band.

## 10. Comparison: one FTF set of FXTFX24-W7-0x-V.3.1 (Q1 – Q3 – Q2)

Source: that notebook's saved Q3D matrix, cell 33 (Q1, Q2, Q3 with a 410 µm coupler
arm). It is re-analysed with this note's method by `simulation/reference/analyse_fxtfx24.py`,
and the matrix is saved as `results/reference/FXTFX24_cell33_cmat_fF.csv`.

| quantity | FXTFX24 set | this design (v4) |
|---|---|---|
| Ec_FX (C_Σ) | 0.927 GHz (20.9 fF) | 1.31 GHz (14.8 fF) |
| Ec_TC (C_Σ) | 0.196 GHz (98.6 fF) | 0.202 GHz (96.0 fF) |
| J f–t (J_1c, J_2c) | 24.8 / 24.7 MHz | 138 / 139 MHz |
| J f–f (J_12) | 1.26 MHz | 12.2 MHz |
| g f–r (β·V_zpf, per n) | 75 / 80 MHz (β = 0.105; f_r 6.109 / 6.475 GHz) | 150 / 149 MHz (β = 0.20) |
| FX charge line | 0.53 fF, β 0.010, T1 limit 3.1 ms | 0.45 fF, β 0.013, T1 limit 2.0 ms |

Notes:
- **Notebook formula for g.** The notebook's own cell-36 formula,
  e·sqrt(π/(2·50·h))·C~⁻¹/2π applied with the readout pad as node, gives 95 MHz on
  this matrix. That is a different convention from β·V_zpf.
- **Cell 36 printouts.** Cell 36's printed numbers (Ec 0.901 GHz, Ec_TC 142 MHz,
  J 194 MHz) come from a different matrix: the last run of the `run_sweep` in
  cell 35. They do not describe this geometry.
- **Why our f–t coupling is ~5.6× stronger:** the flat paddles sit 20 µm from the FX
  pads and the FX capacitance is smaller.
