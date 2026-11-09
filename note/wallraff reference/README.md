# Wallraff group readout circuit (resonator + Purcell filter) — reference

Source: PhD thesis (ETH Zurich, Wallraff group), section 2.5 "Qubit Readout", pp. 119–127
(screenshots in this folder). Their circuit is the same topology as Design6
(feed line → "hugger" C_pf → λ/4 Purcell filter → "hammer" C_pr → λ/4 readout resonator).

| File | Content |
|---|---|
| `p119_fig2.29_readout_circuit_model.png` | Fig. 2.29: lumped-element model of resonator + Purcell filter + feed line; micrograph with all labels |
| `p120_eq2.119-2.120_couplers_lengths.png` | g of qubit–resonator; definition of C_pr, C_prg (hammer), C_pf, C_pff, C_pfg (hugger); length definitions L_rr, L_rrg, L_pf, L_pfg, L_rre, L_pfe; L_in = L_out ≈ 0.33 nH (3 wire bonds) |
| `p121_eq2.121_Cin_typical_J_kappa.png` | C_in ≈ 50 fF input capacitor, |S21| of a series C; typical J/2π ≈ 20 MHz, κ_p/2π ≈ 50 MHz (the "50 kHz" in the text is a typo), ω_p ≈ ω_r − Λ |
| `p122_eq2.122-2.125_ABCD_readout_resonator.png` | ABCD of the qubit end of the resonator, ABCD → Y conversion, adding C_qbrg and C_AB |
| `p123_eq2.126-2.130_Zrq1_Zrq2.png` | Y → ABCD, Z_rq,1 (t-junction → qubit), Z_rq,2 (t-junction → shorted end) |
| `p124_eq2.131-2.134_coupler_Purcell_filter.png` | hammer coupler Z_pr; Purcell filter ABCD; Z_rc seen by the feed line |
| `p125_eq2.135-2.137_full_ABCD_fit_model.png` | full feed-line ABCD (L_in, C_in, feed-line sections, shunt 1/Z_rc, L_out), S21 = 2/(A+B/Z0+CZ0+D); input–output fit model eq. 2.137 |
| `p126_table2.4_design_parameters.png` | Table 2.4: capacitances, lengths, simulated / measured ω_p, ω_r, J, κ_p for qubits 32–36 |
| `p127_fig2.30_feedline_S21.png` | Fig. 2.30: simulated vs measured feed-line transmission (5 channels, 6.5–7.6 GHz) |

## Circuit model (eqs. 2.122–2.136), as admittances

Every "C to ground" is a shunt; C_AB = air-bridge capacitance to ground (value not given in these pages).
Lengths are electrical lengths βL in the ABCD matrices.

* Resonator, qubit side: qubit node = C_qr in series with C_qg to ground, plus C_qrg + 2 C_AB at the
  resonator end, then line L_rr to the hammer t-junction (+ C_AB) → Z_rq,1.
* Resonator, ground side: line L_rrg, 2 C_AB, line L_rre, short → Z_rq,2.
* Hammer coupler: resonator node (1/Z_rq,1 + 1/Z_rq,2 + C_AB + C_prg) — series C_pr — filter side (C_prg + C_AB) → Z_pr.
* Purcell filter (from the feed line): C_AB + C_pff, series C_pf, C_pfg + C_AB, line L_pf, C_AB, shunt 1/Z_pr,
  C_AB, line L_pfg, C_AB, line L_pfe, short → Z_rc.
* Feed line: L_in, C_in, line L_fdl1, shunt 1/Z_rc, line L_fdl2, …, L_out; S21 = 2 / (A + B/Z0 + C Z0 + D).

## Fit model (eq. 2.137, Heinsoo 2018 / Swiadek 2024)

|S21(ω)| = | B ( 1 − A e^{iφ} · 2iκ_p(ω_r − ω) / [4J² + (2i(ω_p − ω) + κ_p) · 2i(ω_r − ω)] ) |

Internal losses are set to zero. A, B, φ describe the background, depth, and asymmetry.
At ω = ω_r, S21 = B: the resonator itself is a transmission *peak*, with the two hybrid dips on either side.

## Table 2.4 (eighth feed line, 43-qubit device)

| | unit | Q36 | Q35 | Q34 | Q33 | Q32 |
|---|---|---|---|---|---|---|
| C_pf | fF | 48.57 | 45.47 | 37.08 | 45.10 | 37.55 |
| C_pff | fF | −5.72 | −5.88 | −6.23 | −5.90 | −6.21 |
| C_pfg | fF | 98.80 | 93.33 | 78.11 | 92.56 | 79.00 |
| C_pr | fF | 2.50 | 2.77 | 2.42 | 2.5 | 2.66 |
| C_prg | fF | 1.95 | 2.39 | 1.83 | 1.95 | 2.21 |
| C_qr | fF | 12.26 | 17.95 | 12.26 | 17.95 | 12.26 |
| C_qrg | fF | 17.61 | 27.03 | 17.61 | 27.04 | 17.61 |
| C_qc | fF | 10.00 | 10.06 | 10.00 | 10.06 | 10.00 |
| C_qcg | fF | 15.40 | 15.51 | 15.40 | 15.51 | 15.40 |
| C_qg | fF | 105.36 | 99.35 | 105.36 | 99.35 | 105.36 |
| C_Σ | fF | 117.62 | 117.30 | 117.62 | 117.30 | 117.62 |
| L_pf | µm | 1200 | 1200 | 1200 | 1200 | 1200 |
| L_pfg | µm | 1920.45 | 2268.05 | 2202.36 | 1844.14 | 2296.28 |
| L_pfe | µm | 2.625 | 2.625 | 2.625 | 2.625 | 2.625 |
| L_rr | µm | 1303.6 | 1300 | 1303.6 | 1300 | 1303.6 |
| L_rrg | µm | 2489.65 | 2749.33 | 2603.73 | 2316.18 | 2720.55 |
| L_rre | µm | 2.625 | 2.625 | 2.625 | 2.625 | 2.625 |
| ω_p^Sim/2π | MHz | 7250 | 6690 | 7050 | 7440 | 6850 |
| ω_p^Meas/2π | MHz | 7258 | 6705 | 7056 | 7451 | 6865 |
| ω_r^Sim/2π | MHz | 7238 | 6671 | 7039 | 7417 | 6839 |
| ω_r^Meas/2π | MHz | 7227 | 6655 | 7018 | 7394 | 6821 |
| J^Sim/2π | MHz | 15.84 | 15.04 | 15.80 | 15.55 | 15.41 |
| J^Meas/2π | MHz | 20.54 | 18.20 | 20.59 | 16.27 | 17.63 |
| κ_p^Sim/2π | MHz | 35.37 | 40.98 | 52.07 | 63.16 | 59.91 |
| κ_p^Meas/2π | MHz | 55.58 | 28.43 | 53.61 | 47.01 | 60.41 |

Other numbers in the text: C_in = 50 fF; L_in = L_out ≈ 0.33 nH; feed-line sections
L_fdl1…7 = 1474.14, 1643.96, 11943.61, 11027.39, 12360.17, 11387.33, 4847.80 µm.

## Key take-aways for Design6

* **Coupling scale:** C_pr ≈ 2.5 fF gives J/2π ≈ 15 MHz in simulation and 16–21 MHz measured.
  The 2.0 fF hammer of Design6 should therefore give J ≈ 12–15 MHz. A 400 MHz splitting in HFSS is **not** explained by C_j.
* **Hugger orientation:** most of the hugger's ground capacitance is on the **filter** side
  (C_pfg ≈ 80–100 fF, C_pff ≈ −6 fF). In Design6 the frame is on the n0 (feed-line) side (110 fF), and only 14 fF is on the filter side.
* **κ_p ≈ 2–4·J:** κ_p/4 < J, so the two hybrid dips with the ω_r peak between them are visible. This is the paper-like shape.
* **Detuning:** ω_p − ω_r ≈ +10 … +23 MHz. The filter sits slightly above the resonator.
* The model, its calibration and the Design6 version are in `Design6/simulation/wallraff_model.py`.
