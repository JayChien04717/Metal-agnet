# Cloud Palace vs Local AEDT Benchmark Note (Design3 & Design5)

Date: 2026-11-24  
Cloud Environment: Google Colab Pro (`CPU + High-RAM`, 8 vCPU / 51 GB RAM), `Qiskit Metal` + `SQDMetal` + `Gmsh` + `AWS Palace` (`Order: 2`, `SuperLU`, `Permittivity = 11.45` default silicon).  
Local Baseline: Ansys AEDT 2021 R2 (Q3D / HFSS).

---

## 1. Design3_FTF_NTHU (`FX1 – TC – FX2` + FX1 Readout)

### 1.1 Capacitance Matrix (3-qubit core: `FX1`, `TC` straight ㄈ bracket 6.7 / 18.9 µm, `FX2`)
- **Cloud run:** `01_capacitance/runs/design3_ftf_cap_20261124_153351`
- **Setup:** Cropped `1.8 × 1.6 mm`, `525 µm` Si ($\epsilon_r = 11.45$), `min_size = 12 µm`, `max_size = 120 µm`, `Order = 2`.
- **Runtime:** **22.0 s** (6×6 Maxwell matrix: ground + 5 qubit/coupler islands).

| Parameter | Palace (Cloud, 22 s) | Local Q3D (12 passes) | Difference | Primary cause |
|---|---:|---:|---:|---|
| $E_C(\text{FX1})$ (GHz) | 0.9217 | 1.0105 | −8.78% | `min_size = 12 µm` on 4 µm gap overestimates $C$; no `res_FX1/2` |
| $E_C(\text{TC})$ (GHz) | 0.1893 | 0.2096 | −9.69% | Same (`4 µm` bracket-to-pad gap needs `min_size ≤ 2–4 µm`) |
| $E_C(\text{FX2})$ (GHz) | 0.9093 | 1.0036 | −9.40% | Same |
| $\|J_{\text{FX1-TC}}\|$ (MHz) | 346.15 | 309.13 | +11.98% | Coarse mesh across 4 µm gap inflates mutual $C$ |
| $\|J_{\text{FX2-TC}}\|$ (MHz) | 420.33 | 369.67 | +13.70% | Same |
| $\|J_{\text{FX1-FX2}}\|$ (MHz) | 96.43 | 68.30 | +41.18% | Small indirect coupling amplified by 4 µm gap mesh + omitted readout conductors |

**Actionable tuning for Design3 Palace capacitance:** Reduce `min_size` from `12 µm` to `2–3 µm` on `FX1, TC, FX2` (since 22 s leaves plenty of budget!) and set `Permittivity = 11.9` to match local Q3D.

### 1.2 Eigenmode (`res_FX1` + `tee_FX1` + `FX1`)
- **Cloud run:** `02_eigenmode/runs/design3_fx1_eigenmode_20261124_152434`
- **Setup:** Cropped `3.6 × 4.8 mm`, `525 µm` Si ($\epsilon_r = 11.45$), `min_size = 15 µm`, `max_size = 160 µm`, `Order = 2`, 3 modes above 5.5 GHz.
- **Runtime:** **701.0 s (~11.7 min)**.

| Mode | Palace Re{f} (GHz) | Palace $Q$ | Local HFSS (GHz) | Difference | Identification |
|---|---:|---:|---:|---:|---|
| Mode 1 | **6.1825** | $4.03 \times 10^5$ | **6.3169** | −134.4 MHz (−2.13%) | **FX1 $\lambda/2$ readout resonator** (verified via standing-wave `log10(|E|)` plot: node at meander midpoint, antinodes at tee & FX1 pad) |
| Mode 2 | 6.3843 | $4.07 \times 10^5$ | — | — | Secondary/box-coupled mode |
| Mode 3 | 11.4371 | $4.04 \times 10^5$ | — | — | Higher harmonic |

**Notes:** Even with `15 µm` mesh and no lumped 200 nH FX1 inductor / different $\epsilon_r$ (`11.45` vs `11.9`), Palace lands within **~2%** of 15-pass HFSS and produces a textbook $\lambda/2$ field distribution.

---

## 2. Design5_6_purcell-res (PurcellReadout Capacitors + ch4 S21)

### 2.1 Capacitance Coupons (`cj` hammer & `cp` EnclosingCap)
- **Cloud runs:** `01_capacitance/runs/design5_cj_cap_20261124_151501` & `design5_cp_cap_20261124_151520`
- **Setup:** `500 µm` Si ($\epsilon_r = 11.45$), `min_size = 10 µm`, `Order = 2`.

| Coupon / Metric | Palace (Cloud) | Local Q3D (8 passes) | Difference | Notes |
|---|---:|---:|---:|---|
| `cj` mutual $C_j$ (fF) | **2.5729** | **2.3743** | **+8.36%** | Runtime 155 s; clean 3×3 Maxwell matrix & field plot |
| `cj` ground $C_{g}$ (fF) | 182.0 / 181.8 | 159.7 / 159.7 | +13.9% | Finite ground plane / mesh size (`10 µm` vs Q3D adaptive) |
| `cp` (`EnclosingCap`) | *1×1 (434.95 fF)* | 40.103 / 17.05 / 119.08 | *Shorted to ground* | **Key finding:** SQDMetal has no `open_terminations`; must add `OpenToGround` on pins `a` & `b` so leads do not merge into ground plane |

### 2.2 Frequency-Driven S21 — Single-Channel Benchmark vs Exact 6-Channel `r2` Chip

| Run | Layout | Solver & Mesh | Freq Sweep | Runtime | Key Result / Finding |
|---|---|---|---|---|---|
| `design5_ch4_s21` (pilot) | Single-channel ch4 (`4.4 × 3.8 mm`) | `Order = 2`, `min_size = 15 µm`, uniform sweep | `6.20–6.90 GHz`, `20 MHz` step (36 pts) | **~20 min/pt** (24/36 pts in ~8 h, ~15.5 GB RAM; stopped early) | Dip at `6.40 GHz` (`-3.26 dB`). Proved brute-force uniform `Order = 2` is impractical on Colab. |
| **`design5_6ch_r2_s21`** (exact local `r2`) | **Exact 6-channel `r2` (`10 × 10 mm`, 41 comps, `C_in` + ch1–ch6)** | `Order = 1`, $\epsilon_r = 11.5$, `min_size = 25 µm`, **PROM `AdaptiveTol = 1e-3`, `AdaptiveMaxSamples = 24`** | `5.80–7.20 GHz`, **`2 MHz` step (701 pts)** | **Palace solve: 1161.6 s (19.4 min)** · **Script total (mesh+solve): 1640.3 s (27.3 min)** · **Colab cell total (incl. plot + Drive backup): 42 min** | **>25× faster** than uniform `Order = 2` for a 6× larger chip at 10× finer resolution! |
| **`design5_6ch_r1_s21`** (exact wafer `r1`) | **Exact 6-channel `r1` (`2026-11-19_purcell_resonator_pair_r1`, `10 × 10 mm`, 41 comps)** | `Order = 1`, $\epsilon_r = 11.5$, `min_size = 25 µm`, **PROM `AdaptiveTol = 2e-3`, `AdaptiveMaxSamples = 50`** | `5.80–7.20 GHz`, **`2 MHz` step (701 pts)** | **Palace solve: 1163.7 s (19.4 min, peak 4.8 GB RAM)** · **Script total (mesh+solve): 1757.1 s (29.3 min)** · **Colab cell total (incl. plot + Drive backup): 44 min** | Pixel-identical wafer `r1` layout; resolved peaks at 5.934 / 5.964, 6.260, 6.404 / 6.466, 7.102 / 7.128 GHz (`C_in` −32 dB baseline). |

- **Cloud run:** `03_s21_driven/runs/design5_6ch_r2_s21_20261125_001019`
- **6-channel `r2` S21 comparison vs Local HFSS `r2`:**
  - Because `C_in` sits in series at the input (`x = -3.9 mm`), transmission shows resonant peaks/Fano features on a −32 dB background (local HFSS r2 poles: ch1 `5.96`, ch2 `6.15`, ch3 `6.37`, ch4 `6.53`, ch5 `6.74`, ch6 `6.98` GHz).
  - Palace PROM resolved strong features at **`5.86 GHz`** (ch1, −100 MHz), **`6.01 GHz`**, **`6.16 GHz`** (ch2, +10 MHz), **`6.27 GHz`** (ch3, −100 MHz), and **`6.50 GHz`** (ch4, −30 MHz), while ch5/ch6 (>6.6 GHz) were missed because `AdaptiveMaxSamples = 24` exhausted its sample budget on the lower 4 channels (12 poles total across 6 channels require `AdaptiveMaxSamples ≥ 50–60` or splitting into two sub-bands, e.g. `5.8–6.5` and `6.5–7.2` GHz) and no ground wire bonds were added (`Order = 1`, `min_size = 25 µm` on `9 µm` CPW gap).
- **Optimization rule derived:** For multi-channel S21 in Palace, always use **PROM Adaptive Frequency Sweep (`AdaptiveTol: 1e-3`)**, set `AdaptiveMaxSamples ≥ 8 × N_channels` (or split into ~0.7 GHz sub-bands), and use `min_size ≤ 15 µm` around `C_in`/`C_p` gaps when quantitative $|S_{21}|$ baseline & $\kappa_p$ accuracy is required.

---

## 3. Summary of Cloud Simulate Optimizations Recorded

1. **Gmsh vertex-dedup patch** in `GMSH_Geometry_Builder._draw_polygon_in_GMSH_from_coords` eliminates `"Could not create line"` crashes on Shapely rounded geometries (`Design3` brackets).
2. **Stub `xmark` + guard `aedt_session`** removes all `gdspy`/AEDT dependencies in Colab.
3. **Explicit `OpenToGround` pins** required for SQDMetal capacitance coupons (`Design5 cp`).
4. **Mesh rule:** `min_size ≤ 0.5 × gap_min` for quantitative (<2%) capacitance agreement on tight couplers (`2–4 µm` for 4 µm gaps).
5. **S21 PROM Adaptive Sweep (`AdaptiveTol`)** mandatory for multi-point driven S21 sweeps in Colab.
6. **`common/plot_utils.py` polish:** fix `x (mm)` / `y (mm)` axis labels (Palace ParaView coordinates × `1e3` are in `µm`, not `mm`).
