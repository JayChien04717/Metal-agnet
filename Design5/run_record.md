# Design5 run record

Structure, functions of each coupling / capacitor, the 6-channel chip and the know-how: see `README.md`.

Design5 is the Design6 W1 variant: one readout channel with a Purcell filter, using the Wallraff-group parameters of
thesis Table 2.4, Q35 (see `SQL HFSS/note/wallraff reference/`).
Targets: ω_r = 6.700 GHz, ω_p = ω_r + 19 MHz, J ≈ 15 MHz, κ_p ≈ 40 MHz (matched-line definition).

## Layout (`pattern design/design5_layout.py`, default = W1)

| | value |
|---|---|
| C_in | enclosing cap, plate 140 µm, next to n0 (feed_c 20 µm) |
| C_p (hugger) | type C: frame (big ground C) toward the **filter**, plate toward n0; plate 165 µm |
| C_j (hammer) | 25 × 65 µm heads, gap 19 µm, pocket margin 16 µm |
| filter | 3.4542 mm, C_p pin → short; tap 1.200 mm from C_p (L_pf) |
| resonator | 4.4070 mm, open end → short; tap 1.300 mm from the open end (L_rr) |
| files | `gds/design5.gds`, `fig/design5_layout.png` |

`D5_VARIANT=P1` builds the old Design6 c6_P1 geometry. `DLF_EXTRA` / `DLR_EXTRA` (mm) lengthen the filter / resonator.

## Method: what is recorded for every run

Run with `simulation/hfss_chan5.py <tag> [f0 f1 count type]`. It runs the HFSS driven modal sweep and then
`analyze5.py`, which appends a section to this file. The run's files go to `simulation/results/<time>_<tag>/`.

1. **Capacitances (C)**
   * Series C_in, C_p, C_j come from the Q3D of the **whole channel in one run** (`q3d_full5.py`: all
     capacitors, lines and pads together; the shorted ends are rendered open so no net merges with ground).
   * The π shunts (electrode → ground) come from the single-capacitor runs. In the full run, a net's ground C
     also contains the line capacitance, which the TL model already includes.
   * Check: the full-run ground C minus L·C' (163.5 fF/mm) equals the shunts of the part runs (filter net 132 fF = 119 + 13).
   * `CAP_SOURCE=part` uses only the part runs.
2. **HFSS S21 and S11**
   * Defaults: Fast sweep (Interpolating misses narrow dips), 0.25 MHz step, basis 2, max ΔS 0.003, 6–16 passes.
   * Recorded: dips, plot, number of passes, last ΔS.
3. **Theory with the Q3D C (no fit):** ABCD model `theory5.py` with the capacitances above and ideal lines (Z0 50.9 Ω, v 1.2016e8 m/s).
4. **Theory fit by C**
   * Free parameters, fitted to the complex S21 and 0.3·S11 over the sweep (differential evolution, then least squares):
     * filter and resonator length corrections;
     * scales of the series C_in, C_p, C_j;
     * scale of v.
   * Soft priors pull the capacitance scales toward 1 (Q3D).
   * Output: fitted C and lengths.
5. **ω_p, ω_r, J, κ_p from a circuit (rows 3 and 4)**, without curve fitting (`modes5.extract`):
   * ω_r = real zero of the branch admittance (S21 = background, the peak between the dips);
   * the two complex poles λ1, λ2 = eigenvalues of [[ω_p + iκ_p/2, J], [J, ω_r]].
   * κ_p is taken **in the channel**: the branch sees C_in on the P1 side, so κ_p is about 1.8× the matched-line value. The matched-line value (Wallraff Table 2.4 definition) is listed under the table.
6. **Direct fit of S21 with eq. 2.137** (Heinsoo 2018, Swiadek 2024), with the background B e^{i(φ_b − 2πfτ)}(1 + b1(f − f_c)):
   * fitted to conj(S21), because HFSS uses e^{+jωt} and eq. 2.137 uses e^{−iωt};
   * parametrised by the two S21 zeros (the lossless nulls fix them sharply);
   * the background is fitted first, then the resonance, then everything together.
7. **Comparison table:** ω_p, ω_r, J, κ_p from rows 3, 4 and 5.

Pipeline test on synthetic data (theory5 with known errors and noise):
* the fit by C recovered all six parameters (dL to 0.2 µm, C scales to 0.1 %);
* the direct fit gave ω_r within 2.4 MHz, J within 1.2 MHz, κ_p within 3 MHz (ω_p within 12 MHz). This is the limit of the 2-mode model.

## Capacitance record

### Q3D, whole channel in one run — q3d5_full (2026-11-09, 5.8 min)

| | part runs (each capacitor alone) | whole channel together |
|---|---|---|
| C_in (P1 feed ↔ node) | 40.57 fF | **41.14 fF** |
| C_p (node ↔ filter) | 45.41 fF | **45.91 fF** |
| C_j (filter ↔ resonator) | 2.59 fF | **3.43 fF** |
| P1 feed ↔ filter | — | 1.04 fF |
| node ↔ resonator | — | 0.03 fF |
| P1 feed ↔ resonator | — | 0.02 fF |

* C_j rises from 2.59 to 3.43 fF (+0.84 fF): the filter and resonator lines add coupling (parallel section, T junctions, stubs).
  That is 32 % more than the hammer alone, but **far too small** to explain the 410 MHz splitting of Design6 c6_P1.
* Ground C beyond the line C (C' = 163.5 fF/mm):
  * filter net 132 fF (part runs: C_p frame 119 + C_j 13);
  * resonator net 21 fF (C_j 13, plus the open end and bends);
  * node net 122 fF;
  * P1 net 207 fF (launch pad).
* Model with these values, at the current lengths:

  | C source | C_j | J | ω_r | ω_p | κ_p in channel | κ_p matched |
  |---|---|---|---|---|---|---|
  | part runs | 2.59 fF | 14.1 MHz | 6.6999 GHz | 6.7297 GHz | 63.6 MHz | 36.6 MHz |
  | whole channel | 3.43 fF | **18.7 MHz** | 6.6930 GHz | 6.7227 GHz | 64.7 MHz | 37.2 MHz |

# Runs
