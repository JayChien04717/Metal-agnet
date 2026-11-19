# Design8: PurcellReadout module test and validation

**Figures:** [fig/](fig/README.md) holds the whole chip ([`chip_fluxonium.png`](fig/chip_fluxonium.png), [`chip_open.png`](fig/chip_open.png), [`chip_transmon.png`](fig/chip_transmon.png)) and one close-up per part, rendered from the current GDS by `code/general/plot_design_figs.py`.

Test chips for the shared Purcell filter + readout resonator module
[`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py), and the Q3D / HFSS runs that check its two
design rules:

1. **Frequency from length.** The filleted length plus the Q3D end loads give f_p and f_r
   ([`purcell_readout_model.py`](../pattern%20design/resonator/purcell_readout_model.py)).
2. **Coupling from the hammer only.** J follows from the Q3D hammer C_j. The layout keeps every other
   filter–resonator path small ([parasitic note](../note/purcell_readout_parasitics.md)).

Status: draft, revision `2026-11-18_r1`. This is a simulation and test geometry. No fabrication GDS, final selection
or release confirmation exists ([design.json](design.json)).

## 1. Layout

One channel on a 5.6 × 4.8 mm chip, 10/6 µm CPW, HFSS library silicon (ε_r 11.9), f_target 6.5 GHz:
feed line → n0 tee → C_p (EnclosingCap, Design5 `ENC_P`) → `PurcellReadout` (400 µm line spacing, taps 1.2 / 1.3 mm,
hammer 25 × 65 µm heads, 30 µm gap) → qubit.

| Build | Resonator end | GDS (simulation geometry, junction-free) | Figures |
|---|---|---|---|
| `open` | OpenToGround; S21 / eigenmode fixture | [gds](gds/Design8_2026-11-18_r1_open_sim.gds) | chip (`Design8_2026-11-18_r1_open_layout.png`, removed 2026-11-19; in git history) · hammer (`Design8_2026-11-18_r1_open_hammer.png`, removed 2026-11-19; in git history) · end (`Design8_2026-11-18_r1_open_qubit.png`, removed 2026-11-19; in git history) |
| `fluxonium` | `MyFluxonium` readout pad | [gds](gds/Design8_2026-11-18_r1_fluxonium_sim.gds) | chip (`Design8_2026-11-18_r1_fluxonium_layout.png`, removed 2026-11-19; in git history) · hammer (`Design8_2026-11-18_r1_fluxonium_hammer.png`, removed 2026-11-19; in git history) · qubit (`Design8_2026-11-18_r1_fluxonium_qubit.png`, removed 2026-11-19; in git history) |
| `transmon` | `MyTmonCoupler` with `readout_line_options.make_rol` | [gds](gds/Design8_2026-11-18_r1_transmon_sim.gds) | chip (`Design8_2026-11-18_r1_transmon_layout.png`, removed 2026-11-19; in git history) · hammer (`Design8_2026-11-18_r1_transmon_hammer.png`, removed 2026-11-19; in git history) · qubit (`Design8_2026-11-18_r1_transmon_qubit.png`, removed 2026-11-19; in git history) |

Build checks ([build_checks_2026-11-18_r1.json](parameters/build_checks_2026-11-18_r1.json)):
- **Pins:** in all three builds, the module pins meet their partners with zero offset, opposite normals and equal widths.
- **Spacing:** no side-by-side CPW pair is closer than 400 µm (skill.md §11).
- **Metal width:** no metal is below 4 µm outside the exempt resonator parts.
- **Lengths:** the drawn filleted lengths equal the targets.

## 2. Validation results (open fixture, revision 2026-11-18_r1)

The lengths come from the frequency model with the Q3D loads below (filter 4.0644 mm, resonator 4.3615 mm, target
f_p = f_r = 6.500 GHz). Nothing was fitted.

**Q3D** ([run 20261118_232405_q3d_module_b3d273d0](results/resonator/20261118_232405_q3d_module_b3d273d0/summary.json),
8 passes, no case reported solver convergence, last-pass change 0.8–4.9 %):
- Hammer C_j: 1.971 fF.
- Head excess: 1.43 fF per side (the two sides give 0.26 / 2.60 fF; mesh noise).
- C_p: 39.37 fF mutual, 14.37 fF plate ground.
- Full-module static mutual: 2.076 fF. The non-hammer part is 0.105 fF, 5.3 % of C_j.

| | Model (Q3D in) | HFSS eigenmode, 6 passes | HFSS driven, 3 passes | HFSS driven, 6 passes |
|---|---|---|---|---|
| f_p (GHz) | 6.5000 | 6.5100 (+10.0 MHz) | 6.4877 (−12.3 MHz) | 6.5082 (+8.2 MHz) |
| f_r (GHz) | 6.5000 | 6.5158 (+15.8 MHz) | 6.4825 (−17.5 MHz) | 6.4906 (−9.4 MHz) |
| J (MHz) | 13.19 (hammer only) | 16.71 (×1.27) | 13.11 (×0.99) | **13.39 (×1.02)** |
| κ_p (MHz) | 26.8 | 47.2 (×1.76) | 25.3 (×0.95) | 25.6 (×0.96) |
| Last ΔF / ΔS | — | 0.157 % | 0.752 | 0.203 |

Runs:
- Eigenmode: [20261118_232801_hfss_eigen_49fd6d7b](results/resonator/20261118_232801_hfss_eigen_49fd6d7b/summary.json).
- Driven, 3 passes: raw [20261118_234933_hfss_driven_0214a698](results/resonator/20261118_234933_hfss_driven_0214a698/run.json),
  analysis [20261119_000312_s21_poles_ce46c76a](results/resonator/20261119_000312_s21_poles_ce46c76a/summary.json).
- Driven, 6 passes: raw [20261119_003221_hfss_driven_73fb5083](results/resonator/20261119_003221_hfss_driven_73fb5083/run.json),
  analysis [20261119_011236_s21_poles_ca7ae9f4](results/resonator/20261119_011236_s21_poles_ca7ae9f4/summary.json).
  The pole fit inside the raw runs is superseded by these re-analyses.
- A 10-pass driven attempt ([20261119_000236_hfss_driven_68062651](results/resonator/20261119_000236_hfss_driven_68062651/run.json))
  ran out of memory at pass 9. Its 8 solved passes had ΔS 0.31–0.93 with no trend, the known skill.md §10 behaviour
  on narrow lossless resonances. The 6-pass budget matches the eigenmode mesh (~525k elements). It is an explicit
  deviation from the 3-pass rule, recorded in run.json.

Reading:
- **J:** the driven J agrees with the Q3D hammer-only closed form within 2 %. So the hammer sets the coupling and
  Q3D is enough to design J.
- **Frequencies:** f_p and f_r are within ±10 MHz (6 passes). The mesh changes them by up to 20 MHz between 3 and
  6 passes, and the eigenmode moves 0.16 % (≈10 MHz) per pass. Frequency agreement is therefore at the level of
  the mesh uncertainty. Trim with `trim_slope` (≈ −1.47 MHz/µm) only after a matching HFSS run.
- **Eigenmode:** it gives 27 % more J and 1.8× κ_p. The κ_p excess is the known eigenmode effect
  (method note §5); the J excess is not explained. Use driven S21 for J and κ_p.
- **Validation status:** unknown (no converged mesh). The qubit-loaded fluxonium / transmon builds were not
  simulated.

## Reproduction

From the repository root, in the `qmetal` environment (`C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe`, AEDT 2021 R2):

```powershell
python -I -B Design8/layout/build.py
python -I -B Design8/simulation/resonator/q3d_module.py
python -I -B Design8/simulation/resonator/hfss_module.py eigen  --q3d Design8/results/resonator/<q3d run>
python -I -B Design8/simulation/resonator/hfss_module.py driven --q3d Design8/results/resonator/<q3d run>
```

Geometry inputs are in [parameters/geometry.json](parameters/geometry.json). `model_loads` there are first estimates.
The HFSS runs replace them with the Q3D values of the named upstream run and recompute the lengths.

## Using the module on a chip

1. Choose `line_spacing` (≥ 400 µm), the taps and the hammer. Run `q3d_module.py` (cases `cj,ref,cp`) for C_j,
   the head excess and the C_p loads.
2. Get the lengths from `purcell_readout_model.length_for_f(f, tap, C_end, s/2, C_j + head_excess)`. For the
   resonator, C_end is the Q3D self-capacitance of the qubit's readout pad, with the qubit present.
3. Get J from `j_hammer(...)` and κ_r from `kappa_readout(J, κ_p, Δ)`. The hammer gap or height for a target C_j comes from `gap_for_cj` /
   `height_for_cj`, a fit to 16/9 µm Q3D data that must be confirmed by the coupon Q3D.
4. Abut the qubit's readout pin to `qubit` and C_p's pin b to `cp`. Put air bridges at `strap_points()`.
