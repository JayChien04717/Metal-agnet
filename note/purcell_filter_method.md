# Transmission Purcell filter (Design4) - calculation methods

Device: input line - interdigital C - lambda/2 CPW (open-open) - interdigital C - output line.

1. Bandwidth: Q_F = f0 / kappa_F. Both ports load the filter equally, so each port has Q_c = 2 Q_F.
2. Coupling capacitor from the lab formula (same as the readout kappa, FXTFX24 cell 42):
   Q_c = pi / (Z0 Zr (Cc w)^2)  ->  Cc = sqrt(pi / (Q_c Z0 Zr)) / w,  Z0 = Zr = 50.9 Ohm, w = 2 pi f0.
   f0 = 6.7 GHz, kappa_F = 200 MHz: Q_F 33.5, Q_c 67, Cc = 101.1 fF.
3. Interdigital capacitance, first estimate (Alley / Bahl), N fingers, overlap l, finger width W, substrate h = 525 um:
   C [pF] = (eps_r + 1) l [(N - 3) A1 + A2],  l in um,
   A1 = 4.409e-6 tanh(0.55 (h/W)^0.45),  A2 = 9.92e-6 tanh(0.52 (h/W)^0.5) pF/um.
   N = 10, W = gap = 6 um, eps_r = 11.45: C = 0.0507 fF/um * l -> l = 199.1 um for 101 fF.
   Replace by Q3D (C_12 of the cap alone) once the layout is confirmed.
4. Filter length: lambda/2 with both ends loaded by Cc: f = v / (2 L) (1 - 2 Cc / (C' L)), C' = 163.4 pF/m,
   v = 1.2016e8 m/s -> L = 7.486 mm. The meander is iterated on the filleted length (as for the readout
   resonators); the metal of the two caps is not counted, so f0 must be checked with HFSS (hfss_scale).
Layout: `Design4/pattern design/purcell_layout.py` (functions alley_cap_fF, coupling_cap, res_length).

# Design6 (Heinsoo / Cai 2024 individual filter) - fitting HFSS S21 with the circuit model

Script: `Design6/simulation/fit_theory.py <tag> [--eig f1,f2]`; model: `theory_abcd.py` (s21_theory, branch_Y).

1. Model: lossless ABCD of P1 - TL - C_in(pi) - TL - n0 - TL - P2 with the shunt branch at n0:
   stub - C_p(pi) - filter (lambda/4, open at C_p, shorted) - tap at D_cj - C_j(pi) - resonator (lambda/4).
   Pi network = Q3D series C + the two electrode-to-ground C (C_ii + C_ij). Z0 = 50.9 Ohm, v = 1.2016e8 m/s.
2. Free parameters: dLf, dLr (effective length corrections, mm), s_in, s_p, s_j (scales on the Q3D series C;
   s_in / s_p also scale the ground C), kv (scale on v). Soft Gaussian priors: s_in, s_p sigma 0.15, s_j 0.4, kv 0.03.
3. Residual: complex S21 and 0.3 x complex S11 in linear scale (not dB, so the lossless -60 dB nulls do not dominate).
   Optional eigenmode constraint: the model eigenfrequencies with the feed cut at n0 (both n0 half lines open,
   0.15 mm each) are the zeros of Im Y(n0) (Foster: B crosses - -> +); residual (f_model - f_HFSS) / 10 MHz.
4. Optimiser: differential_evolution (global) then least_squares.
5. Derived numbers: kappa_p = FWHM of the |S21|^2 notch with C_j -> 0; J = half the minimum splitting of the
   two eigen-cut modes when dLr is scanned through the anticrossing.
   Mode identification: shift dLr / dLf by 20 um and see which mode moves.
Caveat: when the resonator is far detuned (Delta >> J) the data constrain J only through Delta shift J^2/Delta
(< 1 MHz), so J then comes from the prior (Q3D C_j). J is measured only once the two dips are resolved near resonance.

## Wallraff-group circuit model and the W1 design (2026-11-09)

The reference pages are in `note/wallraff reference/` (thesis sec. 2.5, Table 2.4). The script is `Design6/simulation/wallraff_model.py`.

* **Circuit:** eqs. 2.122–2.134, written as admittances. Branch admittance Y_rc at the feed-line T; single channel S21 = 2 / (2 + Z0 Y_rc).
* **ω_p, ω_r, J, κ_p without a curve fit** (eq. 2.137 = 2-mode model):
  * ω_r = real zero of Y_rc (S21 = B, the peak between the two dips);
  * l1, l2 = complex poles of S21 (zeros of 2 + Z0 Y_rc, complex secant) = eigenvalues of [[ω_p + iκ_p/2, J], [J, ω_r]];
  * ω_p = Re(l1 + l2) − ω_r, κ_p = 2|Im(l1 + l2)|, J² = (ω_p + i·Im(l1+l2)) ω_r − l1 l2. Im J² ≈ 0 is the consistency check.
  * Sign: the Y model uses e^{+jωt}, so the poles have Im > 0. Using ω_p − iκ/2 gives Im J² ~ 10⁵: wrong.
* **Reproducing Table 2.4:** one global v = 1.147e8 m/s (their CPW, not given), C_AB = 0.
  * J: 14.0–14.9 vs 15.0–15.8 MHz in the table, within 10 %.
  * ω_p, ω_r: within about 30–90 MHz.
  * κ_p: 27–51 vs 35–63 MHz. The table values come from the full feed line (C_in, standing waves), which the single-channel model leaves out.
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

## Design5 6-channel chip (2026-11-09)

`Design5/README.md` covers the folder structure, the function of each coupling and capacitor, the chip (6.3–7.3 GHz,
200 MHz spacing, one feed line, one C_in), and the know-how list.

Chip lengths: `solve_chip6.py`, using the same capacitors and solving ω_r = f_k, ω_p = f_k + 19 MHz.
Routing: deterministic serpentines in `chip6_layout.snake_points`. The length is exact; each 90° corner subtracts r(2 − π/2).
