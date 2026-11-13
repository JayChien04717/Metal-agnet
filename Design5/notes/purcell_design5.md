# Design5 Purcell-filter notes (6-channel chip)

Design-specific record. General methods: `../../note/purcell_filter_method.md`. Parameters: `../parameters/Design5_<date>.json`.

## W1 design (first Design5 channel, 2026-11-09, superseded by the Design4-style chip)
* **W1 (Design6 CPW, Z0 50.9 Ω, v 1.2016e8 m/s, Q35 capacitances, no qubit / open resonator end, L_pf 1.2 mm, L_rr 1.3 mm):**
  * Lengths for ω_r = 6.700 GHz, ω_p = ω_r + 19 MHz: L_pfg 2426 µm, L_rrg 3155 µm → filter 3.629 mm, resonator 4.458 mm.
  * Result: J 16.9 MHz, κ_p 37.6 MHz; dips at 6689.9 / 6728.4 MHz with a 0 dB peak at 6700 (paper-like).
* **Geometry differences from c6_P1:**
  * the hugger frame must be on the filter side (C_pfg ≈ 93 fF, C_pff ≈ −6 fF; c6_P1 has 14 fF on the filter side and 110 fF on the n0 side);
  * C_pf 45.5 fF (c6_P1 40.6);
  * hammer C_pr 2.77 fF in a small pocket (C_prg 2.4 fF; c6_P1 2.0 fF with 13 fF to ground);
  * taps at 1.2 / 1.3 mm (c6_P1 1.02 mm).

## Design5 pipeline: per-run record (2026-11-09)

Design5 = W1. Every HFSS run is recorded in `Design5/run_record.md` by `Design5/simulation/analyze5.py`, which
`hfss_chan5.py` calls. Each record holds:
* the capacitances, from the whole-channel Q3D (`q3d_full5.py`, all capacitors together);
* HFSS S21 and S11;
* theory with the Q3D C;
* the theory fit by C (lengths, C scales, v);
* the direct eq. 2.137 fit of S21;
* a table comparing ω_p, ω_r, J and κ_p across the three methods.

Notes:
* Whole-channel Q3D: render the shorted ends open (no net merges with ground). Use only the mutual C. The net ground C includes the line C; subtract L/(Z0 v) as a check.
* For the direct eq. 2.137 fit to the HFSS / ABCD complex S21 (e^{+jωt}), fit conj(S21). Otherwise κ runs negative and the fit fails.
* Parametrise the fit by the two S21 zeros z1, z2 (ω_p = z1 + z2 − ω_r, J² = ω_p ω_r − z1 z2), so the sharp lossless nulls are fixed from the start.
* For comparison with a channel S21 fit, κ_p from the circuit must use the real environment at n0 (C_in on one side): `theory5.n0_env_Y`. This is about 1.8× the matched-line κ_p.

## First 6-channel chip, 6.3-7.3 GHz (2026-11-09, superseded)

`Design5/README.md` covers the folder structure, the function of each coupling and capacitor, the chip (6.3–7.3 GHz,
200 MHz spacing, one feed line, one C_in), and the know-how list.

Chip lengths: `solve_chip6.py`, using the same capacitors and solving ω_r = f_k, ω_p = f_k + 19 MHz.
Routing: deterministic serpentines in `chip6_layout.snake_points`. The length is exact; each 90° corner subtracts r(2 − π/2).

## Current chip: six Design4 (Cai 2024) channels, 6.0-7.0 GHz (2026-11-12 / 13)

* Layout `pattern design/chip6_d4.py`: 10 x 10 mm, pads at x = +-4.5 mm, feed line at y = 0, C_in at x = -3.9 mm,
  3 channels up / 3 mirrored down, staggered (feed order ch2, ch4, ch6, ch1, ch3, ch5). Every CPW pair >= 400 um;
  closest channel-to-channel metal 0.879 mm.
* Lengths: `purcell_6500.tune_lengths` at f x 6.08/5.975 (Design4 HFSS-vs-model correction): filter 4.516 -> 3.823 mm,
  resonator 4.832 -> 4.136 mm for 6.0 -> 7.0 GHz; model J 12.7-15.9 MHz, kappa_p 34-49 MHz.
* Full-chip eigenmode (`simulation/chip6_d4_full.py eig`, 12 modes, auto bonds + 54 straps): pass 1 (1.2 M elements)
  gave the pairs 5.907/5.942, 6.104/6.157, 6.317/6.347, 6.514/6.538, 6.694/6.735, 6.882/6.932 GHz: every channel
  60-100 MHz below target. Pass 2 ran out of memory (24 GB) twice, also with 30 um mesh on the lines.
  Next: shorten by f_meas / f_target, tune per channel pair, keep a thinned full-chip run only as a check.
* Full-chip Q3D without bonds / straps (`chip6_d4_full.py q3d`): the g_wb net was not fully connected with bonds.
* Fab GDS `pattern design/chip6_d4_fab.py` -> `gds/Design5_<date>.gds`: lab cheesing, corner test junctions, die marks.
