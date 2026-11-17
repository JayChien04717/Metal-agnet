# Design6 parameter overview (collected 2026-11-17)

This note collects the recorded Design6 Hamiltonian, readout, coupler and flux parameters. The complete r9 parameter snapshot awaits user confirmation (see `wafer_candidate.release_confirmation` in [design.json](../design.json)).

Status labels:
- **target:** design goal.
- **input:** user or fab parameter.
- **derived:** fit or calculation.
- **result:** EM or quantum-model output.

Geometry tag:
- **rev-A:** Historical results for the 2026-11-14 core geometry with flux-line variant A. B2 and r9 require their own runs.
- **r9:** parent wafer revision 2026-11-16_wafer_r9.
- **r10:** current wafer revision 2026-11-17_wafer_r10, with the four trimmed loaded readout tails.

E/h is in GHz and f is not ω. Coupling g is given either in charge units (H term g n_i n_j) or as an exchange (σ⁺σ⁻) value, as stated in each row.

## 1. Qubits and coupler

| Parameter | Q1 | Q2 | TC | Role / status | Source |
|---|---|---|---|---|---|
| f01 at sweet spot (target) | 4.80 GHz | 4.40 GHz | 6.30 GHz | target | [report §2](report_2026-11-15.md) |
| f01 at sweet spot (r9, bare) | 4.81332 GHz | 4.41090 GHz | 6.45108 GHz | derived with drawn L_J | [fixed-L_J model](../results/qubit/20261117_191400_readout_qubit_r9_3cfc15e5/summary.json) |
| f01 at idle (r9, dressed) | 4.80700 GHz | 4.40607 GHz | 6.45919 GHz | r9 Q3D Hamiltonian with drawn L_J | [CZ results](../results/qubit/20261117_181834_cz_r9_20f32abc/cz_summary.json) |
| α at sweet spot (r9, bare) | −204.96 MHz | −207.02 MHz | −212.33 MHz | derived with drawn L_J | fixed-L_J model |
| f01 at idle, dressed (model) | 4.794 GHz | 4.396 GHz | 6.307 GHz | historical rev-A, scqubits on the HFSS fit | [cz_summary.json](../results/qubit/2026-11-14/cz_hfss_v8/cz_summary.json) |
| f01 at Φ0/2 (bare) | 2.68 GHz | 2.45 GHz | 0.75 GHz (d = 0) | derived, historical rev-A | report §4 |
| anharmonicity α (idle, dressed) | −202 MHz | −205 MHz | — | result, historical rev-A | cz_summary.json |
| E_C (HFSS core fit) | 0.185 GHz | 0.185 GHz | 0.188 GHz | derived, historical rev-A (fit RMS 8.2 MHz, 6 passes) | [fit.json](../results/qubit/2026-11-14/hfss_core_v8/fit.json) |
| E_C (Q3D v8) | 0.197 GHz | 0.196 GHz | 0.209 GHz | derived, historical rev-A (Δ 3.1171% at pass 8) | [v8 summary](../results/qubit/2026-11-14/q3d_v8/summary.json) |
| E_J,Σ / E_J,a : E_J,b | 16.84 GHz, 2 : 1 (d = 1/3) | 14.25 GHz, 2 : 1 (d = 1/3) | 28.08 GHz, 1 : 1 | historical rev-A model input (retargeted using the HFSS-fitted E_C) | report §4 |
| E_J/E_C at max / min | 91 / 30 | 77 / 26 | 150 / 0 | derived, historical rev-A | report §4 |
| L_J fab (a / b) | 14.56 / 29.11 nH | 17.21 / 34.41 nH | 11.64 nH each | input (fab target, in `LJ` of the layout) | [design6_layout.py](../layout/design6_layout.py) |
| L_J drawn in the EM runs | 15.89 / 31.77 nH | 18.78 / 37.56 nH | 12.85 nH each | input to historical rev-A EM runs | report §4 |
| I_c (a / b) / R_n (a / b) | 22.6 / 11.3 nA, 12.5 / 25.0 kΩ | 19.1 / 9.6 nA, 14.8 / 29.6 kΩ | 28.3 nA, 10.0 kΩ each | derived (Ambegaokar–Baratoff, Δ_Al 180 µeV) | report §4 |

The EJ values were retargeted from the HFSS E_C on 2026-11-14. Whether this retargeting was user-authorised as final is **unknown**: B2 was confirmed as a candidate only.

## 2. Couplings

| Parameter | Value | Definition | Role / status | Source |
|---|---|---|---|---|
| g Q1–TC / Q2–TC (HFSS fit) | 47.5 / 47.2 MHz | charge units, 4e²C⁻¹_ij/h | derived, historical rev-A | fit.json |
| g Q1–TC / Q2–TC (Q3D v8) | 55.1 / 55.1 MHz | charge units | derived, historical rev-A | v8 summary |
| g Q1–TC, exchange at idle | ≈ 90 MHz | exchange | derived, historical rev-A | report results summary |
| C mutual Q1–TC / Q2–TC | 3.25 / 3.25 fF | Q3D | result, historical rev-A | v8 summary |
| g Q1–Q2 and direct mutual C (Q3D v8) | 2.70 MHz and 0.0556 fF | charge coefficient and Maxwell mutual | historical rev-A; g was used as a prior in the HFSS fit | v8 summary |
| Static ZZ at idle (all at sweet spots) | −5.0 kHz | E11 − E10 − E01 + E00 | result, historical rev-A | cz_summary.json |
| Readout g qubit–resonator (Q3D v8) | 102.4 / 99.3 MHz (Q1 / Q2) | charge coefficient `g/2pi` from the saved claw extraction | derived, historical rev-A | v8 summary |
| Readout g qubit–resonator (earlier report) | 129 / 119 MHz (Q1 / Q2) | historical report value | historical rev-A, calculation behind these values unknown | report §3 |
| XY line → qubit C / T1 limit | 0.085 / 0.084 fF, 287 / 348 µs | Q3D line capacitance and derived T1 | historical rev-A | v8 summary |
| XY crosstalk (other qubit / coupler) | 7.5% / 19% | direct mutual-C ratios | derived, historical rev-A | report §3 |

### Core Q3D on r9

The r9 core fixture completed 20 passes in AEDT 2021 R2 at 5 GHz, Medium order, with a 0.1% target. AEDT reported Δ = **0.27546% at pass 20**, so the run did not converge. The lowest Δ was 0.14899% at pass 19. From pass 19 to 20, the largest changes in E_C and g_n were 0.203% and 0.170%, respectively. All 20 matrices are saved in fF. Pass 20 matches the final export, with seven signal nets and one ground net. The ground-referenced signal matrix is symmetric, has positive diagonals and negative off-diagonals, and is positive definite. [Run record](../results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/run.json) · [Summary and per-pass metrics](../results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/summary_out.json) · [Convergence](../results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/raw/convergence.csv) · [Raw matrices](../results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/raw/)

| Parameter | v8, historical rev-A, pass 8 | r9, pass 8 | r9, pass 20 |
|---|---:|---:|---:|
| Solver Δ | 3.1171% | 3.1385% | 0.27546% |
| E_C/h (Q1 / TC / Q2), GHz | 0.19654 / 0.20898 / 0.19639 | 0.19630 / 0.20879 / 0.19670 | **0.18607 / 0.19708 / 0.18607** |
| C_Σ (Q1 / TC / Q2), fF | 98.55 / 92.69 / 98.63 | 98.67 / 92.77 / 98.47 | **104.10 / 98.28 / 104.10** |
| g_n (Q1–TC / Q2–TC / Q1–Q2), MHz | 55.13 / 55.13 / 2.702 | 55.24 / 55.59 / 2.688 | **52.57 / 52.57 / 2.524** |
| Mutual C (Q1–TC / Q2–TC / Q1–Q2), fF | 3.2518 / 3.2546 / 0.05555 | 3.2654 / 3.2795 / 0.05340 | 3.4730 / 3.4734 / 0.05405 |
| Readout g, charge coefficient (Q1 / Q2), MHz | 102.36 / 99.27 | 102.47 / 99.78 | 106.81 / 103.76 |
| Readout g, legacy transmon convention (Q1 / Q2), MHz | 81.21 / 82.08 | 81.24 / 82.57 | 82.54 / 83.60 |
| XY to own qubit mutual C (Q1 / Q2), fF | 0.08522 / 0.08439 | 0.08829 / 0.09064 | 0.08441 / 0.08572 |
| XY T1 limit (Q1 / Q2), µs | 287.3 / 347.7 | 268.0 / 300.9 | 310.0 / 356.5 |
| XY1 drive crosstalk (Q2 / TC) | 8.45% / 22.29% | 7.68% / 24.68% | 8.38% / 24.65% |
| XY2 drive crosstalk (Q1 / TC) | 8.37% / 22.39% | 7.69% / 24.38% | 7.89% / 24.12% |
| XY1 direct C ratio (Q2 / TC) | 7.49% / 18.65% | 6.89% / 21.12% | 7.59% / 21.18% |
| XY2 direct C ratio (Q1 / TC) | 7.43% / 18.73% | 6.94% / 20.85% | 7.15% / 20.64% |

The v8 raw matrix reproduces its saved summary under the same [analysis code](../simulation/qubit/q3d_core.py). At equal pass count, r9 and v8 differ by at most 0.16% in E_C and 0.84% in g_n. Continuing r9 from pass 8 to 20 lowers E_C by 5.2–5.6% and Q1–TC/Q2–TC g_n by 4.8–5.4%. The larger differences between r9 pass 20 and v8 pass 8 therefore cannot be attributed to the geometry change alone. [v8 raw matrix](../results/qubit/2026-11-14/q3d_v8/cmat_fF.csv) · [v8 convergence](../results/qubit/2026-11-14/q3d_v8/convergence.csv) · [r9 pass 8](../results/qubit/20261117_161936_q3d_core_r9_3b00f1f5/raw/pass_8_fF.csv)

The r9 fixture uses the wafer XY routes and grounded B2 coupler flux geometry. It renders the three islands, both claws with their 40 µm stubs, and both XY routes with open ends. Launch pads, separate qubit flux traces, and the Purcell channels are outside the selected render. E_C and g_n follow from the inverse island capacitance block with the other conductors grounded. Readout g uses a 50 Ω λ/4 approximation at the 7.00/6.80 GHz target frequencies. The legacy transmon-convention g and 50 Ω XY T1 use the helper's harmonic estimate at the target qubit frequencies. Drive crosstalk uses the ratio of inverse-matrix voltage-division coefficients. Direct C ratios use the Maxwell mutuals.

The user selected pass 20 as the current r9 design input on 2026-11-17 and accepted approximately 3% design accuracy. Adoption records which numbers to use. AEDT's Δ measures adaptive change and the absolute extraction error remains unknown. Final E_J authorization and full-device validation remain pending. [Design index](../design.json)

### XY termination candidates

The adopted r9 extraction gives a 50 Ω drive-port T1 limit of 310/357 µs and 10 MHz Rabi power of −84.10/−83.87 dBm on chip (Q1/Q2). Adding the drive-port decay rates to independent 100 µs intrinsic decay gives 75.6/78.1 µs total T1.

Using the same Q3D capacitance but exact bare sweet-spot f01 and |⟨0|n|1⟩| at the drawn L_J gives 320/370 µs and −83.94/−83.70 dBm (Q1/Q2). [Fixed-L_J model](../results/qubit/20261117_191400_readout_qubit_r9_3cfc15e5/summary.json)

The local Q3D fixture retains the r9 qubits, coupler, claws and terminal bends and cuts the XY lines at x = ±1.7 mm. Shortening both terminal segments increases the offset from the horizontal qubit-arm etch. All points use pass 3 and a 3% solver stop target. [Comparison run](../results/qubit/20261117_183814_xy_gap_compare_r9_249bd090/summary.json)

| Offset | Own mutual C, Q1 / Q2 | XY T1, Q1 / Q2 | 10 MHz Rabi, Q1 / Q2 | TC drive ratio, XY1 / XY2 | Other qubit, XY1 / XY2 |
|---|---|---|---|---|---|
| 250 µm, local baseline | 0.0963 / 0.1046 fF | 0.220 / 0.222 ms | −85.59 / −85.93 dBm | 24.0% / 23.4% | 4.84% / 4.13% |
| 350 µm | 0.0518 / 0.0509 fF | 0.760 / 0.934 ms | −80.20 / −79.69 dBm | 30.5% / 33.9% | 6.73% / 7.24% |
| 425 µm | 0.0395 / 0.0478 fF | 1.307 / 1.062 ms | −77.85 / −79.13 dBm | 35.2% / 32.5% | 7.64% / 6.62% |

At 425 µm, the absolute TC voltage-division coefficient drops by 39%/36% from the local baseline. Its ratio to the own-qubit drive rises because the own coupling drops faster. Local baseline C is 14%/22% higher than the adopted full-route result. Q1/Q2 own mutual C changed by 42%/38% from pass 2 to 3 at 425 µm despite AEDT Δ = 0.928%. Pass-to-pass C drift at 350 µm is unknown because the first intermediate matrix triggered a zero-denominator error after the final matrix was exported. The 425 µm geometry is a screening candidate and its full-route T1 and crosstalk remain unknown.

## 3. Readout (Purcell filter + λ/4 resonator)

| Parameter | Q1 | Q2 | Role / status | Source |
|---|---|---|---|---|
| f_r = f_p target (loaded, with claw) | 7.00 GHz | 6.80 GHz | target | `READOUT` in design6_layout.py |
| L_filter / L_res (fab, r9 parent) | 3.82263 / 4.03611 mm | 3.94469 / 4.15878 mm | historical geometry input | [design.json](../design.json) `parent_readout_lengths_mm` |
| L_filter / L_res (fab, r10) | 3.880590 / 4.053848 mm | 3.993377 / 4.137843 mm | selected geometry input | [design.json](../design.json) `readout_lengths_mm` |
| `CLAW_DL` (resonator shortening for the claw) | 0.100 mm (≈ 17 fF) | same | input | design6_layout.py |
| Taps: filter / resonator | 1.20 / 1.30 mm | same | input | layout |
| CPW (readout) | 16 / 9 µm, Z0 49.9 Ω, v 1.199e8 m/s | same | input / derived (conformal) | layout |
| Hammer (r9) | 200 µm centres, 25 × 37.5 µm head, gap 40 µm, ground 16 µm | same | input | design.json `readout_hammer_um` |
| C_p (coupon / full readout Q3D) | 41.25 / 42.17 fF | 41.25 / 42.17 fF | result; coupon not converged (0.14 %); full Q3D not converged | [readout note](readout_s21_20261117.md) |
| Hammer C_j (coupon, deembedded) | 0.66 fF | 0.66 fF | result, converged | readout note |
| Filter–res static mutual (full Q3D) | 2.19 fF | 2.23 fF | result, not converged (stable) | readout note |
| **HFSS r9 readout-only** (+0.1 mm L_res, no claw): bare f_p / f_r | 7.114 / 7.121 GHz | 6.904 / 6.929 GHz | result, wide grid, ΔS 0.0091, validation unknown | readout note (`cf91f657`, `af4e9eaf`) |
| HFSS r9: J / κ_p | 16.6 / 38.2 MHz | 16.1 / 33.9 MHz | result, same run | readout note |
| Q3D-corrected circuit model (variant G): f_p / f_r / J / κ_p | 7.111 / 7.136 GHz, 18.2 / 36.6 MHz | 6.905 / 6.932 GHz, 17.6 / 32.0 MHz | derived (no HFSS fitting), validation unknown | readout note (`07c856b2`) |
| Length candidate ΔL_filter / ΔL_res | +65.96 / +80.74 µm | +66.69 / +83.06 µm | variant G candidate | [loaded comparison](readout_s21_20261117.md#length-candidate-and-loaded-s21) |
| Loaded candidate HFSS f_p / f_r / J / κ_p | 6.98585 / 6.89244 GHz, 15.41 / 35.72 MHz | 6.76919 / 6.63358 GHz, 13.95 / 32.03 MHz | three adaptive passes, ΔS = 0.061393 | loaded comparison |
| Loaded candidate κ_r/2π and model 2χ | 0.873 MHz, −2.114 MHz | 0.324 MHz, −1.668 MHz | resonator-like HFSS pole FWHM and isolated bare-resonator model shift | loaded comparison |
| r10 loaded f_p / f_r / J / κ_p | 6.99175 / 6.98892 GHz, 15.96 / 36.15 MHz | 6.80975 / 6.82477 GHz, 15.25 / 31.66 MHz | recorded loaded fixture, three adaptive passes, ΔS = 0.070565 | [before and after trim](readout_s21_20261117.md#loaded-hfss-tail-trim-candidate) |
| χ / readout target | ≈ 0.7 MHz, κ_r ≈ 2χ | ≈ 0.7 MHz | target | report §2 |
| rev-A HFSS + pyEPR (drawn L_J, qubits 4.44 / 4.02 GHz): 2χ, κ_eff | 1.50 MHz, 1.64 MHz | 1.14 MHz, 1.70 MHz | result, rev-A, historical | report §6 |
| Purcell T1 limit | ≫ 1 ms (below eigenmode resolution) | same | result, historical rev-A | report §6 |

Open readout points:
- **Filters above target:** the r9 filters sit about 100 MHz above target (HFSS f_p 7.114 / 6.904 vs 7.00 / 6.80 GHz).
- **Loaded candidate detuning:** the tail trim brings f_r − f_p to −2.83/+15.02 MHz from −93.4/−135.6 MHz (Q1/Q2). The wafer-length loaded operating point remains unknown.

## 4. Two-qubit gate (CZ)

| Scheme | Parameters | Result | Role / status | Source |
|---|---|---|---|---|
| Q1 + TC pulsed together, \|11⟩↔\|20⟩ | Q1 → 4.575 GHz, TC → 5.058 GHz, 8 ns ramps | **80 ns: F 0.99962, leakage 1.1e-4, φ 180.06°** (coherent) | result, rev-A, historical | cz_summary.json |
| Coupler-only, historical baseline targets | idle ZZ −5 kHz, ZZ_on −0.84 MHz at f_c = f_Q1 + 0.3 GHz | 700 ns: F 0.99996 | historical rev-A exploration | [CZ coupler note](cz_coupler_flux_20261117.md) |
| Coupler-only, historical straddling proposal | Q1 / Q2 / TC 4.55 / 4.40 / 6.80 GHz, g_qc 61 MHz, g12 4.0 MHz | idle ZZ +5.4 kHz, 160 ns: F 0.99949, leakage 1.5e-4 | historical rev-A exploration | CZ coupler note |
| Q1 + TC, r9 drawn L_J | idle ZZ −16.80 kHz, f_Q1 / f_TC at gate 4.5832 / 5.0479 GHz | 60 ns: F 0.999254, leakage 1.16e-4 | coherent model result | [r9 comparison](cz_coupler_flux_20261117.md#r9-hamiltonian) |
| Coupler-only, r9 drawn L_J | idle ZZ −16.80 kHz, gate f_TC 5.2081 GHz | 700 ns: F 0.999986, leakage 1.18e-5 | coherent model result | r9 comparison |
| Coupler-only, r9 E_C, hypothetical straddling and native g_n | idle ZZ +17.53 kHz | 200 ns: F 0.973681, phase 203.87° | coherent model result at drive bound | r9 comparison |
| Coupler-only, r9 E_C, hypothetical straddling and g_qc / g12 = 61 / 4 MHz | ZZ zero at f_TC 6.6324 GHz | 160 ns: F 0.998295, leakage 1.42e-3 | coherent model result with hypothetical inputs | r9 comparison |

## 5. Flux lines (B2, current layout)

| Parameter | Q1 | TC | Q2 | Role / status | Source |
|---|---|---|---|---|---|
| Geometry | L type, sep 9 µm, bar 45 µm, fbl_x 12 µm | same | same (mirrored) | input (B2 confirmed by the user 2026-11-15) | README §3 |
| Mutual inductance M | 1.26 pH | 1.34 pH | 1.26 pH | result, B2 (raw-run mapping incomplete) | report header |
| DC crosstalk max / rms | 5.9 % / 4.0 % | | | result, B2 | report header |
| Pulse distortion (on chip, variant A) | ≤ 0.016 % ring-up within 30 ns | | | result, historical rev-A | report §7 |

## 6. What is still missing for a complete, current parameter set

- **Core Hamiltonian on r9:** the adopted Q3D values are recorded above. The r9 HFSS core fit remains pending.
- **Final EJ authorisation:** the user decides whether to keep the retargeted L_J (14.56 / 29.11, 17.21 / 34.41, 11.64 nH) or change the qubit targets. The straddling CZ proposal, for example, needs Q1 f_max 4.55 GHz.
- **Loaded readout:** the r10 fixture retains 2.83/15.02 MHz absolute detuning (Q1/Q2). The solver missed its ΔS target. Hybrid-mode χ remains unknown.
- **CZ scheme choice:**
  - Q1 + TC pulse, native coupler-only, or the shorter coupler-only proposal with modified g_qc and g12.
  - The current r9 pulse comparison is recorded above. Dissipation and sensitivity checks remain open.
