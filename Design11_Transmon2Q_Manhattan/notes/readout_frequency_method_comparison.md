# Readout Frequency Method Comparison (Q3D Model vs HFSS Eigenmode vs HFSS Driven S21)

Date: 2026-11-28. Scope: Design11 PurcellReadout channels Q1 / Q2. Current channel geometry is `r17` (unchanged through `r20`):

| Channel | `L_filter` | `L_res` | `head_gap` | `C_p` plate | Hammer spacing |
|---|---:|---:|---:|---:|---:|
| Q1 | 4.1577 mm | 4.4060 mm | 45.2 µm | 267.0 µm | 200 µm |
| Q2 | 4.5024 mm | 4.5370 mm | 44.6 µm | 123.0 µm | 200 µm |

All frequencies are in GHz and all differences in MHz. "Δ vs target" is method minus target. Pairwise differences are row minus column, e.g. "Eig − Q3D".

## 1. Current geometry (`r17` channels, `r20` build)

### Source runs

| Method | Run(s) | Settings | Geometry match |
|---|---|---|---|
| **Q3D-calibrated circuit model** | `20261126_155318_q3d_length_solve_11ac23ec` | Q3D line constants (`20261126_133418_q3d_cpw_line_cpw106_a7600d2d`) + `r16` in-situ loads (`20261126_143615_q3d_channel_insitu_r16_full_6f7f7bec`) + `r16` model check (`20261126_144132_q3d_model_check_fa91743e`) | `r17` channels, solved directly. The loads use the old core with `cross_length` 143.4 / 143.5 µm, so they are `needs_revalidation` for `r18`+ (§3, item 6). |
| **HFSS eigenmode** | Solve `20261128_132318_hfss_eigen_r17_6p_5e4a52bd` → re-analysis `20261128_164443_eigen_from_sd_r17_pass3_5e03bacb` | `hfss_readout_eigen.py`, explicit 10 µm mesh, min f 5.9 GHz, 6 modes. Launched at the 6-pass default. On 2026-11-28 the user said 3 passes are enough, so it was stopped during pass 4. The pass-3 (LastAdaptive) eigenvalues were read from the AEDT `*.sd` text files with `simulation/resonator/eigen_from_sd.py`. AEDT was not reopened. | yes (current `build.py`) |
| **HFSS driven S21** | Solve `20261128_125428_purcell_s21_driven_r17_4p_272cdb6a` → fit `20261128_132341_fit_s21_r17_4p_b889cde1` | `hfss_readout_s21.py --mode driven --passes 4`, explicit 10 µm mesh, Fast sweep 5.75–6.75 GHz in 1 MHz steps. Complex fit with `--window Q2=6.04,6.21 --window Q1=6.29,6.45`. | yes |

### Q1 (target 6.35 GHz)

| Method | `f_p` | `f_r` | Δ`f_p` vs target | Δ`f_r` vs target | `f_r − f_p` | `J` (MHz) | `κ_p` (MHz) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Target | 6.3500 | 6.3500 | 0 | 0 | 0 | 4.743 | 30.0 |
| Q3D model | 6.3500 | 6.3500 | 0.0 | 0.0 | 0.0 | 4.743 | (input 30) |
| Eigenmode, pass 3 | **6.3844** | **6.3650** | +34.4 | +15.0 | −19.4 | 4.90 | 31.1 |
| Driven S21, fitted pole | 6.3999 | 6.3448 | +49.9 | −5.2 | −55.1 | 5.26 | 30.9 |
| Driven S21, `\|S21\|` filter dip | 6.373 | — | +23 | — | — | — | — |

| Pair | Δ`f_p` | Δ`f_r` |
|---|---:|---:|
| Eig − Q3D | +34.4 | +15.0 |
| Driven − Q3D | +49.9 (dip +23) | −5.2 |
| Driven − Eig | +15.5 (dip −11) | −20.2 |

### Q2 (target 6.15 GHz)

| Method | `f_p` | `f_r` | Δ`f_p` vs target | Δ`f_r` vs target | `f_r − f_p` | `J` (MHz) | `κ_p` (MHz) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Target | 6.1500 | 6.1500 | 0 | 0 | 0 | 4.743 | 30.0 |
| Q3D model | 6.1500 | 6.1500 | 0.0 | 0.0 | 0.0 | 4.743 | (input 30) |
| Eigenmode, pass 3 | **6.1560** | **6.1649** | +6.0 | +14.9 | +8.9 | 5.09 | 30.3 |
| Driven S21, fitted pole | 6.1542 | 6.0736 | +4.2 | −76.4 | −80.6 | 4.27 | 30.0 |
| Driven S21, `\|S21\|` filter dip | 6.164 | — | +14 | — | — | — | — |

| Pair | Δ`f_p` | Δ`f_r` |
|---|---:|---:|
| Eig − Q3D | +6.0 | +14.9 |
| Driven − Q3D | +4.2 (dip +14) | −76.4 |
| Driven − Eig | −1.8 (dip +8) | −91.3 |

### Eigenmode details (observed; not evidence of convergence)

**Per-pass record** (`EMConvData`). Max ΔF: pass 2 47.0 %, pass 3 3.11 %. Solved elements: 519,138 → 622,802 → 747,277. `IsConverged = false`. Pass 4 was interrupted, and its data was not used.

**Raw pass-3 modes** (f, κ = 2f''):

| Mode | κ / Q | Use |
|---|---|---|
| 6.1060 GHz | κ ≈ 4.96 GHz, Q ≈ 1.2 | port / box mode, dropped |
| 6.1551 GHz | κ = 27.9 MHz | Q2 filter-like |
| 6.1658 GHz | κ = 2.48 MHz | Q2 resonator-like |
| 6.3642 GHz | κ = 1.14 MHz | Q1 resonator-like |
| 6.3852 GHz | κ = 29.9 MHz | Q1 filter-like |
| 11.83 GHz | — | out of band |

**Bare frequencies.** Each hybridised pair is inverted with `purcell_resonator_pair_sim.invert`, giving the bare `f_p`, `f_r`, `J` and `κ_p` above.

**Pass-to-pass movement** of the bare values:

| Quantity | Pass 1 → 2 | Pass 2 → 3 |
|---|---:|---:|
| Q1 `f_p` | +52.6 MHz | +7.9 MHz |
| Q1 `f_r` | +71.5 MHz | +33.0 MHz |
| Q2 `f_p` | +48.2 MHz | +2.4 MHz |
| Q2 `f_r` | +68.5 MHz | +17.9 MHz |

The filter frequencies are nearly settled. The resonator frequencies were still rising by 18–33 MHz per pass, so treat eigenmode `f_r` as uncertain by at least that amount, with a likely upward bias still to come.

### Driven details (observed; not evidence of convergence)
- **Convergence:** ΔS per pass 0.326 → 0.097 → 0.061 (4 passes, 915,629 elements, 26.6 min). This is above the 0.02 target, so the run is not converged.
- **Q1 fit:** relative residual 0.30 %. The Fano amplitude is η = 2.05 and the phase φ = −1.06 rad. η > 1 cannot occur in an ideal hanger; it absorbs the `C_in` standing wave. As a result, the Q1 pole and the |S21| dip disagree by 27 MHz.
- **Q2 fit:** relative residual 2.1 %. The Q2 resonator feature at 6.072–6.073 GHz is sampled by 1–2 points, so `J` is uncertain to about ±1 MHz.
- **Superseded fit:** `20261128_132234_fit_s21_r17_4p_ce59b978` used the default Q2 window (6.08–6.21 GHz), which excludes the resonator. Its Q2 values are invalid.

## 2. Historical, `r16` channels (not valid for `r17`)

`r16` differs from `r17` only in channel lengths and `head_gap`: Q1 `L_filter` +25.1 µm, `L_res` +4.4 µm and `head_gap` 38.9 µm; Q2 `L_filter` −4.0 µm, `L_res` −16.1 µm and `head_gap` 38.8 µm.

| Method | Source run | Q1 `f_p` / `f_r` | Q1 `J` / `κ_p` | Q2 `f_p` / `f_r` | Q2 `J` / `κ_p` |
|---|---|---:|---:|---:|---:|
| Q3D-calibrated model | `20261126_144132_q3d_model_check_fa91743e` | 6.3148 / 6.3432 | 5.53 / — | 6.1542 / 6.1703 | 5.51 / — |
| Eigenmode, pass 2 of 3 (stopped by user during pass 3; modes read approximately from the AEDT plot) | `20261126_144828_hfss_eigen_r16_nowb_3p_52819a4d` | coupled modes 6.344, 6.358 | — | coupled modes 6.143, 6.160 | — |
| Driven 3 passes (ΔS 0.126) | `20261126_140926_…_r16_nowb_3p_7d4253b4` → fit `20261126_143332_fit_s21_r16_nowb_3p_fe66850b` | 6.3518 / 6.3540 | 5.66 / 28.6 | 6.1498 / 6.0902 | 5.38 / 30.6 |
| Driven 2 passes (ΔS 0.384) | `20261126_134423_…_r16_nowb_2p_cf1afc8a` → fit `20261126_140020_fit_s21_r16_nowb_2p_b5cdb28d` | 6.3418 / 6.3480 | 5.76 / 28.5 | 6.1394 / 6.1238 | 5.44 / 30.4 |
| Driven 2 passes with wire bonds (parent status `failed` after the S-parameter export) | `20261126_130307_…_r16_s200_2p_fast_clean_ec655610` → fit `20261126_131703_fit_s21_r16_s200_2p_fast_clean_7b1abf5b` | 6.3558 / 6.3595 | 5.60 / 29.2 | 6.1317 / 6.0809 | 5.08 / 28.8 |

**Driven pass-to-pass scatter at `r16`:** Q2 `f_r` spans 6.081–6.124 GHz (43 MHz). Q1 `f_r` spans 6.348–6.360 GHz (12 MHz).

## 3. Interpretation

1. **`κ_p` agrees between eigenmode and driven within 1 MHz** (31.1 / 30.9 MHz for Q1, 30.3 / 30.0 MHz for Q2), and both are on the 30 MHz target.
   - **J:** eigenmode gives 4.90 / 5.09 MHz and driven 5.26 / 4.27 MHz, against the 4.74 MHz target. The two HFSS methods differ by 0.4–0.8 MHz. Driven Q2 `J` rests on 1–2 sweep points.
   - **Use:** driven S21 is usable for κ_p. J is uncertain to about ±0.5 MHz.
2. **Driven Q2 `f_r` is unreliable.** It is 91 MHz below eigenmode and 76 MHz below the Q3D model, and the `r16` driven runs show the same −80 MHz offset.
   - **Cause:** this matches the repo lesson that driven S21 resonator frequencies are unreliable, because the adaptive mesh follows the feed and filter modes.
   - **Consequence:** driven Q1 `f_r` is also 20 MHz below eigenmode. Do not take `f_p` / `f_r` from driven S21.
3. **Q2 filter is consistent across all three methods.** `f_p` (Q3D 6.1500, eigenmode 6.1560, driven pole 6.1542) agrees within 6 MHz.
4. **Q1 filter is high in both HFSS methods:** eigenmode +34 MHz, driven pole +50 MHz, driven dip +23 MHz. At `r16`, driven was +37 MHz from the model.
   - **Likely cause:** the Q3D model systematically under-predicts the Q1 filter frequency. A candidate is the filter end load with the large 267 µm `C_p` plate, which is 70.6 fF on Q1 against 38.7 fF on Q2 in the `r16` model check.
   - **Status:** the eigenmode filter frequencies moved only 2–8 MHz from pass 2 to pass 3, so this offset is probably not a mesh artifact. It is unverified beyond 3 passes.
5. **Resonator `f_r`: eigenmode is +15 MHz from the model on both channels.** The pass 2 → 3 movement (+18 / +33 MHz) is larger than this offset, so the eigenmode `f_r` cannot yet confirm or refute the model to better than about ±20–30 MHz.
6. **Geometry caveat.** The Q3D model loads come from `r16` with the old core (`cross_length` 143.4 / 143.5 µm, before the T flux lines). Both HFSS runs use the `r20` build (`cross_length` 150.4 / 150.5 µm). Part of the +15 MHz `f_r` offset may come from the changed claw end load. It has not been quantified.
7. **Summary.**
   - **Q2:** the three methods agree within about 15 MHz on `f_p` and on Q3D / eigenmode `f_r`. The −76 to −91 MHz driven `f_r` is an artifact.
   - **Q1:** `f_p` is 34 MHz above target in eigenmode (pass 3) and `f_r` is 15 MHz above target. That is outside the ±7.5 MHz detuning tolerance (`f_r − f_p` = −19 MHz in eigenmode).
   - **Trim:** any length trim should be decided from a higher-pass eigenmode run (resonator frequencies still moving), not from this 3-pass result.

## 4. Not done / limitations
- **Eigenmode stopped at 3 passes** (user override of the 6-pass default, 2026-11-28). Resonator modes were still moving 18–33 MHz per pass, and max ΔF was 3.11 % at pass 3.
- **Driven not converged:** ΔS 0.061 after 4 passes. Pass counts are execution settings, not convergence evidence.
- **Other eigenmode runs:**
  - `20261128_164008_hfss_eigen_r17_3p_user_ccc7aee2` was stopped by the user (`stopped_by_user`) and has no results.
  - A later 3-pass restart was not executed.
  - The pass-3 data used here comes from the cancelled 6-pass run. Its AEDT results folder `simulation/aedt/D11.aedtresults/E_5e4a52bd_hfss.results` was left untouched.
- **`design.json` adopted runs are unchanged.** `adopted_runs.readout_s21` still points to the `r15` driven run. The `r17` runs above are recorded only in this note and in their run folders.

## 5. Fast single-channel HFSS eigenmode (`r17` channels, `r20` build, 2026-11-28)

To get 6-pass eigenmode frequencies under 1 hour per run, `simulation/resonator/hfss_readout_eigen.py` was extended with opt-in fast-tuning flags (preserving existing full-chip defaults):
- `--channel Q1|Q2`: render only that channel (filter, resonator, its qubit/claw, and local feed section; Q1 keeps `feed_cin` + `cin` + `feed_in` + `n0_Q1` + `feed_mid`; Q2 uses `feed_mid` + `n0_Q2` + `feed_out` without `C_in`).
- `--local-mesh`: apply explicit 10 µm length mesh only to filter/resonator CPW traces, their gaps (non-model sheets), coupling hammer, `C_p`, and qubit readout claw; rest of chip adapts freely.
- `--no-wire-bonds`: disable auto wire bonds on etched paths while tuning (keeps the 6 hammer PEC straps per channel).
- `--n-modes 4 --min-freq <target - 150 MHz>` (6.20 GHz for Q1, 6.00 GHz for Q2).
- `--pass-by-pass --early-stop-mhz 3.0`: read per-pass eigenfrequencies, Q, and solved elements directly from `*.sd` / `*.EMConvData` after each pass and save `pass_history.json` immediately; stop early if every channel mode moves < 3 MHz between passes.

### Source runs

| Channel | Run | Settings | Solver time | Pass-6 elements | Early stopped? |
|---|---|---|---:|---:|---|
| **Q1** | `20261128_173735_hfss_eigen_r17_q1_fast2_773b81a3` | `--channel Q1 --local-mesh --no-wire-bonds --n-modes 4 --min-freq 6.20 --passes 6 --pass-by-pass --early-stop-mhz 3.0` | **26.82 min** | 474,877 | No (pass 5→6 max mode Δ = 11.4 MHz) |
| **Q2** | `20261128_180600_hfss_eigen_r17_q2_fast_719c96a8` | `--channel Q2 --local-mesh --no-wire-bonds --n-modes 4 --min-freq 6.00 --passes 6 --pass-by-pass --early-stop-mhz 3.0` | **23.88 min** | 444,278 | No (pass 5→6 max mode Δ = 6.5 MHz) |

*(An earlier Q1 attempt `20261128_170045_hfss_eigen_r17_q1_fast_0ea7ac6b` finished pass 1 in 2.8 min but hung on `project.Save()` inside the pass loop; superseded by `773b81a3`.)*

### Q1 comparison (target 6.3500 GHz)

| Method | `f_p` (GHz) | `f_r` (GHz) | Δ`f_p` vs target | Δ`f_r` vs target | `f_r − f_p` | `J` (MHz) | `κ_p` (MHz) | Elements / Time |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Target | 6.3500 | 6.3500 | 0 | 0 | 0 | 4.743 | 30.0 | — |
| Q3D model (`11ac23ec`) | 6.3500 | 6.3500 | 0.0 | 0.0 | 0.0 | 4.743 | (input 30) | — |
| **Fast eigenmode Q1, pass 6 (`773b81a3`)** | **6.3868** | **6.3858** | **+36.8** | **+35.8** | **−1.0** | 5.16 | 23.7 | 474,877 / 26.82 min |
| Previous 3-pass full eigenmode (`5e03bacb`) | 6.3844 | 6.3650 | +34.4 | +15.0 | −19.4 | 4.90 | 31.1 | 747,277 / stopped in pass 4 |

**Q1 per-pass movement (`773b81a3`):**

| Pass | Time (min) | Elements | Bare `f_p` (GHz) | Bare `f_r` (GHz) | Δ`f_p` (MHz) | Δ`f_r` (MHz) | Max mode Δ (MHz) | `J` (MHz) | `κ_p` (MHz) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2.62 | 190,870 | 6.31400 | 6.23168 | — | — | — | 8.28 | 1.16 |
| 2 | 2.46 | 229,026 | 6.35407 | 6.31052 | +40.1 | +78.8 | 78.6 | 7.17 | 27.23 |
| 3 | 2.99 | 274,824 | 6.36250 | 6.34379 | +8.4 | +33.3 | 32.8 | 5.74 | 14.02 |
| 4 | 3.93 | 329,777 | 6.36971 | 6.36400 | +7.2 | +20.2 | 19.8 | 5.25 | 18.06 |
| 5 | 4.97 | 395,732 | 6.37909 | 6.37558 | +9.4 | +11.6 | 11.8 | 5.21 | 20.09 |
| 6 | 8.40 | 474,877 | 6.38683 | 6.38582 | +7.7 | +10.2 | 11.4 | 5.16 | 23.70 |

### Q2 comparison (target 6.1500 GHz)

| Method | `f_p` (GHz) | `f_r` (GHz) | Δ`f_p` vs target | Δ`f_r` vs target | `f_r − f_p` | `J` (MHz) | `κ_p` (MHz) | Elements / Time |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Target | 6.1500 | 6.1500 | 0 | 0 | 0 | 4.743 | 30.0 | — |
| Q3D model (`11ac23ec`) | 6.1500 | 6.1500 | 0.0 | 0.0 | 0.0 | 4.743 | (input 30) | — |
| **Fast eigenmode Q2, pass 6 (`719c96a8`)** | **6.1784** | **6.1821** | **+28.4** | **+32.1** | **+3.7** | 5.04 | 11.9* | 444,278 / 23.88 min |
| Previous 3-pass full eigenmode (`5e03bacb`) | 6.1560 | 6.1649 | +6.0 | +14.9 | +8.9 | 5.09 | 30.3 | 747,277 / stopped in pass 4 |

*\*Q2 single-channel omits `C_in`, so its standing-wave enhancement (~2×) is absent and eigenmode `κ_p` is ~12 MHz rather than ~30 MHz; use driven S21 for `J` and `κ_p`, and eigenmode only for `f_p` / `f_r`.*

**Q2 per-pass movement (`719c96a8`):**

| Pass | Time (min) | Elements | Bare `f_p` (GHz) | Bare `f_r` (GHz) | Δ`f_p` (MHz) | Δ`f_r` (MHz) | Max mode Δ (MHz) | `J` (MHz) | `κ_p` (MHz) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2.39 | 178,582 | 6.10466 | 6.08070 | — | — | — | 8.23 | 0.27 |
| 2 | 2.18 | 214,290 | 6.15823 | 6.13977 | +53.6 | +59.1 | 59.4 | 6.81 | 1.08 |
| 3 | 2.80 | 257,132 | 6.16783 | 6.15617 | +9.6 | +16.4 | 16.6 | 5.25 | 2.27 |
| 4 | 3.42 | 308,550 | 6.16976 | 6.16677 | +1.9 | +10.6 | 9.0 | 5.03 | 4.32 |
| 5 | 4.33 | 370,236 | 6.17319 | 6.17467 | +3.4 | +7.9 | 6.1 | 5.02 | 7.24 |
| 6 | 7.45 | 444,278 | 6.17843 | 6.18211 | +5.2 | +7.4 | 6.5 | 5.04 | 11.90 |

### Proposed length trims (not applied; awaiting user confirmation)

Using `df/dL ≈ −f/L` (shorted-arm trim only; hammer spacing stays 200 µm):

| Channel | Line | Current `L` (mm) | `df/dL ≈ −f/L` (MHz/µm) | Pass-6 Δ`f` vs target (MHz) | Proposed Δ`L` (µm) | Proposed `L` (mm) |
|---|---|---:|---:|---:|---:|---:|
| Q1 | `L_filter` | 4.1577 | −1.527 | +36.8 | **+24.1** | 4.1818 |
| Q1 | `L_res` | 4.4060 | −1.441 | +35.8 | **+24.9** | 4.4309 |
| Q2 | `L_filter` | 4.5024 | −1.366 | +28.4 | **+20.8** | 4.5232 |
| Q2 | `L_res` | 4.5370 | −1.356 | +32.1 | **+23.7** | 4.5607 |

### Convergence status and caveats
- **Pass counts are execution settings, not convergence evidence.** Neither run reached the 3 MHz early-stop threshold within 6 passes: at pass 5 → 6, Q1 bare modes still moved +7.7 / +10.2 MHz (max mode Δ = 11.4 MHz, HFSS max ΔF = 6.35 %) and Q2 bare modes moved +5.2 / +7.4 MHz (max mode Δ = 6.5 MHz, HFSS max ΔF = 5.10 %).
- Both channels show nearly zero filter–resonator detuning by pass 6 (`f_r − f_p` = −1.0 MHz on Q1, +3.7 MHz on Q2, well within ±7.5 MHz), confirming that the large −19.4 MHz detuning seen in the 3-pass full eigenmode was mostly the slower pass-to-pass convergence of the resonator relative to the filter.
- Because frequencies are still drifting upward by ~5–10 MHz at pass 6, the proposed trims above represent a lower-bound correction (~21–25 µm longer); `build.py` and `design.json` are unchanged until confirmed.
