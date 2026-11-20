# Design4: PurcellReadout hammer-stub test (153.5 → 75 µm)

**Question.** Can the hammer stub of the `PurcellReadout` module be halved, and does the Q3D circuit model still
predict the HFSS driven S21 when it is?

Status: test geometry, revision `2026-11-20_r1` ([design.json](design.json)). No fabrication GDS, final selection or
release confirmation. The former Design4 (Cai 2024 side-by-side pair) was deleted on 2026-11-20; it is in git history
up to `fb1130e`.

## 1. Geometry

One channel, open fixture (no qubit), on a 5.6 × 4.8 mm chip: feed line (two 50 Ω ports) → n0 tee → C_p (Cai T-type
`EnclosingCap`, plate 182 µm, gap 12 µm) → `PurcellReadout` → `OpenToGround`. Same cross-section and hammer as the
Design5 / Design6 channels: 16/9 µm CPW, ε_r 11.5, heads 25 × 65 µm, head gap 27 µm, taps 1.2 / 1.3 mm,
f_target 7.0 GHz ([geometry.json](parameters/geometry.json)).

| Variant | Metal stub (main-line edge → head outer edge) | `line_spacing` |
|---|---|---|
| `ref` | 153.5 µm (Design5 / Design6 as drawn) | 400 µm |
| `half` | 75 µm | 243 µm |

Simulation GDS and figures: [gds/](gds/), [fig/](fig/) (`build.py`, estimate lengths). Build checks
([build_checks_2026-11-20_r1.json](parameters/build_checks_2026-11-20_r1.json)): pins meet with zero offset, no metal
below 4 µm, drawn filleted lengths equal the targets. The HFSS geometries use the lengths from their own Q3D loads
(saved in each driven run).

## 2. Results

**Q3D** (8 passes; no case reported solver convergence; last-pass change 0.8–5.6 %):

| | C_j (fF) | Head excess (fF) | C_p mutual / plate ground (fF) | Static non-hammer mutual | Run |
|---|---|---|---|---|---|
| `ref` | 2.298 | 0.98 | 40.10 / 17.05 | 0.142 fF (6.2 % of C_j) | [q3d_ref](results/resonator/20261120_123447_q3d_ref_d6860330/summary.json) |
| `half` | 2.411 | 0.90 | 40.10 / 17.05 | 0.243 fF (10.1 % of C_j) | [q3d_half](results/resonator/20261120_123817_q3d_half_f2637ce6/summary.json) |

**HFSS driven against the Q3D circuit model** ([compare run](results/resonator/20261120_132804_compare_0ba94282/summary.json),
[figure](results/resonator/20261120_132804_compare_0ba94282/model_vs_driven.png)). Both S21 traces go through the same
pole extraction; nothing is fitted. HFSS: 6 adaptive passes, one Fast sweep 6.5–7.5 GHz at 1 MHz.

| | J HFSS / model (MHz) | J ratio | f_p, f_r HFSS − model (MHz) | κ_p HFSS / model (MHz) | Last ΔS | Side-by-side J bound |
|---|---|---|---|---|---|---|
| `ref` | 18.73 / 16.30 | ×1.15 | +24.0, +30.4 | 34.6 / 34.3 | 0.216 | ≤ 0.77 MHz |
| `half` | 21.19 / 17.32 | ×1.22 | +26.3, +38.7 | 35.4 / 34.8 | 0.111 | ≤ 1.64 MHz |

What this shows:
- **The halved stub works as a coupler.** The S21 keeps the same two-dip shape. J rises 13 % in HFSS
  (18.7 → 21.2 MHz) and 6 % in the model (16.3 → 17.3 MHz; C_j +5 %).
- **The model underestimates J more at 75 µm.** The HFSS-minus-model gap grows from 2.4 MHz (`ref`) to 3.9 MHz
  (`half`). The extra 1.4 MHz is of the size of the side-by-side T-region coupling at 243 µm (`parallel_j` bound
  rises by 0.9 MHz; static non-hammer mutual 6 % → 10 %), which the hammer-only model leaves out. The 2.4 MHz already
  present at 400 µm is not explained here.
- **κ_p agrees** within 2 % in both variants.
- **Frequencies:** HFSS sits 24–39 MHz (0.3–0.6 %) above the model in both variants. At 6 passes the mesh is not
  converged (ΔS 0.22 / 0.11), and Design8 moved by about 20 MHz between 3 and 6 passes, so this offset is not yet a
  length correction.

Limitations: unconverged Q3D and HFSS meshes (above); 1 MHz sweep step (pole-fit RMS 0.010–0.012 on both traces);
open fixture without a qubit; PEC wire straps, not air bridges; no eigenmode cross-check.

Runs not used: `20261120_123100_q3d_ref_fa8deb77` (failed: module-case net names matched twice; fixed),
`20261120_124233_hfss_driven_ref_fe54b209` (failed before solving: Windows path length in the code snapshot),
`20261120_130425_drv_ref_41a5b91b` (5001-point sweep refused by AEDT). The `ref` HFSS data come from
`20261120_124310_drv_ref_344da20b`: its adaptive solve and 1 MHz Wide sweep completed, then the separate fine sweep
was refused by AEDT, so the run is marked failed but `raw/S_wide.csv` is complete.

## 3. Reproduction

From the repository root, qmetal Python (`C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe`, AEDT 2021 R2):

```bash
python -I -B Design4/layout/build.py
python -I -B Design4/simulation/resonator/q3d_module.py --variant ref
python -I -B Design4/simulation/resonator/hfss_driven.py --variant ref --q3d <q3d_ref run> --step 0.001
python -I -B Design4/simulation/resonator/compare.py <drv_ref run> <drv_half run>
```

Repeat the Q3D and HFSS steps with `--variant half`.
