# ZZ and readout dispersive shift of the FTF set (Design1)

Script: `Design1/simulation/zz_analysis.py [EJ] [EL] [--extras]` (numpy / scipy only, ~1 min).
Log: `Design1/simulation/results/runs/zz_EJ4_EL0.6.log`, map: `Design1/fig/zz/zz_map_J12_J1c.png`.

## Model

Ec and J come from the latest Q3D run (`results/latest/fxarray_fl15_summary.json`):
Ec = 1.322 / 0.202 / 1.326 GHz (FX1 / TC / FX2), J_1c / J_2c / J_12 = -140.0 / +141.9 / -16.5 MHz
(J = 4 e^2 C~^-1[i,j] / h, coefficient of n_i n_j with n the Cooper-pair number operator).

    H = sum_i H_i + J_1c n_1 n_c + J_2c n_2 n_c + J_12 n_1 n_2
    FX:  4 Ec n^2 + EL phi^2 / 2 - EJ cos(phi - pi)      (HO basis 160, 8 levels kept)
    TC:  4 Ec n^2 - EJc cos(phi)                         (charge basis, 6 levels kept; EJc set by f_c)

ZZ = E_101 - E_100 - E_001 + E_000 of the dressed states (largest overlap with the bare product
state, coupler in |0>). ZZ is invariant under the sign / phase gauge of the n operators; only the
sign of J_1c J_2c J_12 matters (here +, which gives the cancellation between the direct and the
coupler-mediated term).

## Result for EJ = 4, EL = 0.6 GHz (half flux)

FX: f01 = 0.507 GHz, f12 = 4.00 GHz, f03 = 7.06 GHz, |n01| = 0.117, |n03| = 0.377. n13 = 0
(parity), so the 6.3 / 6.5 GHz resonators do not hit f13.

| | |
|---|---|
| ZZ (f_c 3 - 10 GHz) | -0.1 ... -0.3 kHz, nearly flat in f_c |
| direct J_12 alone (J_1c = J_2c = 0) | -1.5 kHz |
| coupler alone (J_12 = 0) | -0.6 ... -0.8 kHz |
| ZZ vs J_12 at J_1c = 140 MHz | 0 MHz: -0.67 kHz, 10 MHz: -0.03, 16.5 (design): -0.21, 30: -2.1, 60: -13 kHz |
| EJ +-10 %, EL +-10 % | ZZ -0.15 ... -0.26 kHz, f01 380 - 670 MHz |
| exchange g_|100>-|001> (f_c 3 - 6 GHz) | 0.055 - 0.058 MHz (direct alone 0.226 MHz, cancelled by the coupler) |
| readout chi, FX1 - R1 (g 161 MHz, 6.3 GHz) | 7.3 MHz (kappa 1.2 MHz) |
| readout chi, FX2 - R2 (g 160 MHz, 6.5 GHz) | 8.5 MHz |
| chi at g = 80 / 60 MHz (R1) | 1.8 / 1.0 MHz |

The device sits in the ZZ valley of the (J_12, J_1c) map (same valley as the device of the paper,
J_12 = 40, J_1c = 236 MHz).
