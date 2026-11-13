# Design4 Purcell-filter notes (single channel)

Design-specific record: tuning history, Q3D / HFSS numbers, decisions. General methods: `../../note/purcell_filter_method.md`.
Current design: section "Design4 = Cai 2024 parameters" at the end. Parameters: `../parameters/Design4_<date>.json`.

## First version: transmission filter with interdigital caps (6.7 GHz, superseded)

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

## Ground bonds on every CPW are required (Design4, 2026-11-11)

Rule: in HFSS (and on the chip) every CPW of the readout channel needs ground bonds / air bridges, including
the short pieces: hammer stubs, T-junction arms, the stub to the hugger, meander segments. Without them the
circuit (capacitor) model cannot match HFSS.

* Cause: stub_f + hammer pocket + stub_r form one etched channel that cuts the ground between the filter and the
  resonator. The return current detours and gives a coupling that is not a capacitance; Q3D (ground = one net)
  cannot see it, so an S21 fit puts it all into C_j. Wallraff's model has C_AB (air bridges) at every
  t-junction for the same reason (thesis Fig. 2.29).
* qiskit-metal `auto_wirebonds` only bonds straight segments > 400 um (`wb_threshold`), so all the short
  pieces were unbonded.
* HFSS eigenmode, 250 um spacing, paper hammer gap 9 um (Q3D C_j 2.40 fF -> capacitor model J 16.2 MHz):

  | ground bonds | J (MHz) | equivalent C_j (fF) | kappa_p (MHz) |
  |---|---|---|---|
  | qiskit auto bonds only | 68.2 | 10.3 | - |
  | + 2 straps on the hammer stubs | 22.8 | 3.39 | 66.9 |
  | 62 straps on every CPW (<= 200 um apart) | 20.4 | 3.04 | 61.2 |

  At 190 um spacing the coupling was 122 MHz splitting even with the hammer electrodes removed (continuous pocket).
* Code: `Design4/pattern design/purcell_6500.py` `ground_straps(design, 'hammer' | 'all' | 'cpw')`,
  `draw_hammer_straps(renderer, design, which=...)` after `render_design` (PEC bond wires named g_wb).
  Run: `Design4/simulation/eigenmode/hfss_hammer_straps.py --straps cpw`. With 62 straps AEDT faults on the second
  design of a session (COM -2147417851): one trim per fresh AEDT. Results: `Design4/simulation/results/2026-11-11/`.
* Still open: kappa_p 61 MHz (HFSS eigenmode) vs 39 MHz (circuit model, Q3D C_p 39.8 fF); the bonds do not change it.
* Calibration for design (HFSS / Q3D): J = 8.5 MHz per fF of Q3D C_j (circuit model 6.7 MHz/fF).

## 400 um spacing, target chi/2pi = 1 MHz (transmon), 2026-11-11

Targets: kappa_r = 2 chi = 2 MHz, kappa_p ~ 50 MHz -> J ~ 4.9 MHz at f_p = f_r.

* Q3D, paper hammer, 400 um spacing (`Design4/simulation/q3d/q3d_hammer_gap.py`, results `q3d400_gap_sweep/`);
  C_p 39.76 fF at every gap:

  | gap (um) | 9 | 20 | 30 | 45 | 60 | 80 |
  |---|---|---|---|---|---|---|
  | C_j (fF) | 2.202 | 1.291 | 0.916 | 0.619 | 0.464 | 0.356 |

* No-hammer eigenmode, 400 um, straps on every CPW (`nohammer400_cpwwb/`): direct filter-resonator coupling
  J_par = 2.14 MHz at the crossing (trim -11 um), 1.25-1.29 MHz at +-40 MHz detuning (per-point 2-mode inversion).
  Bonding every CPW does not remove it (1.2 MHz without the extra straps). Likely the parallel runs near the taps.
* The 3-trim scan fit (fixed f_p, linear f_r) fails for weak coupling (gave 6.8 MHz): use the per-point inversion
  (`point_J` in hfss_hammer_straps.py). For the 250 um hammer case both agree (20.4 MHz at the crossing).
* If J_par adds in phase: hammer must give ~2.8 MHz -> C_j 0.41 fF (circuit model 6.73 MHz/fF) -> gap ~69 um;
  0.33 fF (HFSS 8.5 MHz/fF) -> ~88 um (extrapolated). The phase of J_par vs the hammer is NOT verified:
  run the hammer at ~70 um, 400 um spacing, near the crossing, before fixing the gap.

### Hammer at 400 um, gap 70 um, straps on every CPW (2026-11-12)

`hfss_hammer_straps.py --straps cpw --spacing 0.4 --gap 70` (66 bonds), results `Design4/simulation/results/2026-11-12/hammer400_gap70_cpwwb/`.

| trim (um) | modes (GHz) | detuning (MHz) | per-point J (MHz) | kappa_p (MHz) |
|---|---|---|---|---|
| -40 | 6.4990 / 6.5247 | +25.3 | 3.05 | 57.8 |
| -11 | 6.4912 / 6.5048 | -13.3 | 3.70 | 58.4 |
| +20 | 6.4431 / 6.4980 | -54.7 | 2.37 | 55.7 |

* Crossing near trim ~ -24 um (f ~ 6.50 GHz); J at the crossing ~ 3.5 MHz (interpolated) -- below target.
* kappa_p here ~ 57 MHz (not 50): the target for kappa_r = 2 MHz becomes J ~ 5.2 MHz.
* Q3D C_j(70 um) ~ 0.40 fF. Neither phase hypothesis matches the 6.7-8.5 MHz/fF calibration:
  adding J_par gives 3.3 MHz/fF (-> gap ~29-33 um for 4.9-5.2 MHz), subtracting gives 13.9 MHz/fF (-> gap ~53-55 um).
  Both say the gap must shrink from 70 um. A second gap (~40 um) fixes the slope.
* Capacitor model fitted to these eigenmodes (`eigenmode/s21_model_vs_eigen.py`: per trim, L_filt, L_res, C_p, C_j
  solved so the two circuit poles of the full channel with C_in equal the two HFSS modes; figure
  `s21_model_vs_eigen.png`): C_p 49.2 / 49.5 / 48.3 fF (Q3D 39.76, +22-24 %), C_j 0.458 / 0.564 / 0.367 fF
  at trim -40 / -11 / +20 um (Q3D 0.403). The kappa_p excess is an effective C_p of ~49 fF. C_j changes with
  detuning, so one capacitor cannot hold the direct (no-hammer) coupling path.

### Driven S21, 400 um / gap 70 / trim -11 um (2026-11-12)

`Design4/simulation/driven/hfss_hammer_driven.py`, results `results/2026-11-12/hammer400_gap70_driven_m11um/`.
* Mesh Link to the solved eigen design (P400g70_cpwwb_m11um_0141_hfss) was silently refused by AEDT (ImportMesh stays
  false, no message). Fallback: driven adaptive on S at 6.5 GHz, same seeds, 8 passes, first order: 474k tets, Delta S
  still 0.21 (not converged). Fast sweep 6.35-6.65 GHz keeps fields; |E| exported at 6.4912 / 6.5048 GHz.
* Filter: S21 fit gives C_p 38.6 fF, kappa_p 36.4 MHz -- matches Q3D (39.76 fF, 38.5 MHz) and the earlier driven
  filter-only 34.9 MHz, NOT the eigenmode kappa_p ~58 MHz. So the kappa_p mismatch sits in the eigenmode kappa, not in C_p.
* Resonator: no feature anywhere in the sweep. |E| maps show the resonator ~100x weaker than the filter and a coarse
  mesh along it (the adaptive on S never refines the undriven resonator). The resonator / J part of this driven run is
  not usable; it needs the eigenmode mesh (Mesh Link working) or a mesh refined on the resonator.

### S21 only (no eigenmode), 400 um / gap 70 / trim -11 um (2026-11-12)

`driven/hfss_hammer_driven.py` (adaptive on S at 6.49 GHz, seeds 10 / 5 um, first order, 8 passes -> 1.10 M tets,
Delta S 0.037; Fast sweep 6.40-6.60 GHz, 0.05 MHz, fields saved), results `results/2026-11-12/hammer400_gap70_s21_m11um/`.
* Filter: capacitor-model fit to S21 gives C_p 39.22 fF, kappa_p 37.6 MHz; Q3D 39.76 fF -> 38.6 MHz. Agreement 1.4 %.
* Resonator: no line in S21. Q3D C_j (0.403 fF) predicts a clear Fano line near 6.480 GHz. An Interpolating sweep
  6.47-6.50 GHz on the same mesh (`driven/interp_sweep.py`, tol 0.2 %) is also smooth (residual from a smooth curve
  4e-5), so it is not a Fast-sweep artifact. Fitted C_j (0.154 fF) only reflects the missing line, not a measurement.

## Wide layout: every CPW pair >= 400 um apart (2026-11-12)

Rule from Jay: every two side-by-side CPW runs at least 400 um apart, centre to centre.
* The old folded tails broke it in 12 places: meander pitch 230 um, clearance to the main line 250 um,
  feed_jf next to the filter meander 250 um.
* New routing: `purcell_6500.serpentine_tail`; `build(wide=True)` is now the default.
  * Filter: up 400 um from jf, then a vertical run at x = -0.4 mm, then x = -0.8 mm.
  * Resonator: down at x = 0.4 mm, then up at x = 0.8 mm.
  * The lowest runs are 650 um above the feed line.
  * Same drawn lengths (4.1913 / 4.5060 mm), same taps and hammer.
* Check: `cpw_spacing_check.py` (old: 12 violations, wide: none). Figure / GDS: `Design4/fig|gds/purcell_wide400_gap45.*`.
* HFSS / Q3D design names start with `W`.

### Q3D capacitance matrix, wide layout, 400 um, paper hammer

Script: `q3d_hammer_gap.py`, run WITHOUT wire bonds, 20 passes (folder deleted 2026-11-12; superseded by the run with bonds below).

| gap (um) | 30 | 45 | 60 | 70 |
|---|---|---|---|---|
| C_j (fF) | 0.904 | 0.604 | 0.456 | 0.393 |

* C_p is 39.83 fF and C_in 39.92 fF at every gap. The old layout gave C_j 0.40 fF at 70 um, so the meanders hardly change C_j.
* Full matrix at gap 45 um (fF; Q3D nets frame_cin = input, frame_cp = feed, head_a = filter, head_b = resonator):

  | | input | feed | filter | resonator | ground |
  |---|---|---|---|---|---|
  | input | 374.12 | -39.92 | -0.100 | -0.011 | -332.23 |
  | feed | | 660.42 | -39.83 | -0.089 | -578.71 |
  | filter | | | 758.10 | -0.604 | -716.67 |
  | resonator | | | | 758.63 | -757.92 |

* Gap for J target 5.2 MHz (kappa_p ~ 57 MHz -> kappa_r 2 MHz):
  * HFSS calibration 8.5 MHz/fF: C_j 0.61 fF -> gap 44.5 um;
  * circuit model 6.73 MHz/fF: C_j 0.77 fF -> gap 35 um.
  * Gap 45 um was chosen for the S21 run.

## Solver settings from now on (Jay, 2026-11-12)

**Q3D 8 passes; eigenmode 6 passes; S21 (driven) 3 passes; HFSS initial length mesh 10 um on every CPW
line and capacitor electrode** so the mesh starts fine and does not need many refinement passes. These are the
script defaults (`q3d_hammer_gap.py`, `hfss_hammer_straps.py`, `hfss_hammer_driven.py`); see skill.md section 12.
(The wide400 Q3D sweep above ran with the older 20-pass / 0.3 % setting.)

## Ground wire bonds: qiskit-metal hfss_wire_bonds (Jay, 2026-11-12)

* Every etched CPW is flagged `hfss_wire_bonds` (`purcell_6500.set_wire_bonds`). The renderer options
  `P.WB_OPTIONS` = wb_threshold 350 um, wb_size 3 (span 66 um), 15 um PEC wire. This puts one bond per
  350 um of every straight segment: 28 bonds (4 on each ~1.5 mm run, 1 on each 400 um serpentine bottom, 7 on the feed).
* The hammer stubs and the tee arms (20-50 um long) never get auto bonds. They get
  `draw_hammer_straps(renderer, d, 'all')`: 9 straps (2 hammer stubs, n0 x 3, jf / jr arms).
* The same bonds are used in Q3D, eigenmode and driven (script defaults: `--straps all` + auto bonds).
  Q3D counts the `g_wb` nets as ground.
* Figure: `Design4/fig/purcell_wide400_gap45_wirebonds.png` (`pattern design/plot_wirebonds.py`).
* All old results (`Design4/simulation/results/*`, the old AEDT project `simulation/aedt/`) were deleted on 2026-11-12.
  Result paths quoted in the older sections above no longer exist; their numbers are kept as history.

## Design4 = Cai 2024 parameters (Jay, 2026-11-12) — CURRENT DESIGN

Reference: Cai et al., Chin. Phys. B 33, 068501 (2024), `note/wallraff reference/Cai_2024_Chinese_Phys._B_33_068501.pdf`
(password-protected; read with `pdftotext -layout`).

| | Cai 2024 | Design4 now |
|---|---|---|
| f_r = f_p | 6.08 GHz (Fig. 1f) | 6.08 GHz (`P.F0`) |
| CPW w / s | 16 / 9 um (least sensitive to etch dL) | 16 / 9 um |
| C_in, C_p | 40 fF, T-type enclosing cap, L = 12 um | T-type enclosing, gap 12 um, plate 182 um |
| C_j | 2 fF, two rectangular plates | plates 25 x 65 um, gap 42 um |
| J, kappa_p | 8.8, 45.5 MHz | (S21 below) |
| Si eps_r | 11.5 | 11.5 (`SI_EPS` default in purcell_6500) |

* Lengths: `tune_lengths()` at 6.08 GHz (circuit model with Design5 Q3D caps) -> filter 4.536 mm, resonator 4.853 mm
  (`P.LENGTHS`; model J 12.7, kappa_p 33.9 MHz). HFSS decides the final trim.
* Layout: wide serpentines (all CPW >= 400 um), 30 qiskit auto bonds + 9 hammer / tee straps;
  `gds/purcell_cai_6080.gds`, `fig/purcell_wide400_gap26_wirebonds.png`.
* Q3D tuning, 8 passes, with bonds (~2 min each), `results/2026-11-12/cai_q3d*`:

  | plate (um) | C_j gap (um) | C_in (fF) | C_p (fF) | C_j (fF) |
  |---|---|---|---|---|
  | 190 | 26 | 41.26 | 41.18 | 2.779 |
  | 185 | 36 | 40.57 | 40.48 | 2.204 |

  C_j ~ 53.8 / gap + 0.71 fF; C_p ~ 0.14 fF per um of plate width -> plate 182 um, gap 42 um.
* Rule: one Q3D per geometry change, no gap sweeps; two points are enough to interpolate.
  | 182 | 42 (FINAL) | 40.09 | 40.03 | 1.956 |

### Final design: Q3D matrix + HFSS S21 (2026-11-12)

Q3D (8 passes, bonds, `results/2026-11-12/cai_final_q3d/gap42um/`), fF:

| | input (C_in frame) | feed (C_p frame) | filter (C_j head a) | resonator (C_j head b) | ground |
|---|---|---|---|---|---|
| input | 372.22 | -40.09 | -0.195 | -0.019 | -330.58 |
| feed | | 650.62 | -40.03 | -0.175 | -568.61 |
| filter | | | 792.74 | -1.956 | -749.71 |
| resonator | | | | 790.16 | -787.53 |

HFSS driven (3 passes, 10 um mesh, 390k tets, Delta S 0.048, 13.7 min; `results/2026-11-12/cai_final_s21/`).
Figure `cai_fit.png` from `driven/cai_fit_plot.py`: HFSS S21 / S11, the ideal model (Q3D C, drawn lengths), and
the model fitted to |S21| + |S11| (L_filt, L_res, C_in, C_p, C_j free):

| | Cai | ideal (Q3D) | fit to HFSS |
|---|---|---|---|
| C_in / C_p / C_j (fF) | 40 / 40 / 2 | 40.09 / 40.03 / 1.956 | 42.90 / 41.09 / 0.669 |
| L_filt / L_res (mm) | | 4.536 / 4.853 | 4.637 / 4.913 |
| f_p / f_r (GHz) | 6.08 / 6.08 | 6.081 / 6.083 | 5.958 / 6.018 |
| kappa_p / J (MHz) | 45.5 / 8.8 | 33.1 / 11.7 | 33.0 / 4.0 |

* Filter: HFSS dip at 5.9475 GHz, 133 MHz below the design. Equivalent to +0.10 mm of filter length
  -> shorten the filter by ~0.10 mm (resonator by ~0.06-0.10 mm) for 6.08 GHz.
* C_p from the fit is within 3 % of Q3D (C_in within 7 %): the feed-filter coupling is OK. kappa_p is 33 MHz with
  C_p = 40 fF; Cai's 45.5 MHz would need C_p ~ 47 fF (kappa_p ~ C_p^2), or their kappa includes the environment.
* Resonator: no line in the HFSS S21 (the same as the 10/6 runs; with the filter 130 MHz below the resonator it is
  far detuned and narrow). The fitted C_j (0.67 fF) and J are NOT a measurement: they only fit the absence of the line.
  J needs f_p = f_r (fix the lengths first) or an eigenmode run.

### 2026-11-12 hammer width 20 um, gap 30 um (user request)
Q3D 8 passes (results/2026-11-12/q3d_w20_gap30): C_j 2.388 fF, C_p 40.03 fF (was w25/g42: C_j 1.956 fF).
Higher C_j -> J rises above the eigenmode 15.2 MHz (Cai target 8.8 MHz). Lengths not yet re-tuned (both ~0.08 mm too long).
Reverted 2026-11-12: hammer back to width 25 / height 65 / gap 42 um (C_j 1.956 fF); the w20/g30 test is for reference only.


## Implementation in Design4 code (rules from skill.md sections 11-13)

* CPW spacing check: `pattern design/cpw_spacing_check.py` -> `report(design, min_um=400)`.
* Wide routing: `purcell_6500.serpentine_tail`; `build(..., wide=True)` is the default (`wide=False` = old 230 um
  folded tails). HFSS / Q3D design names start with `W` for the wide layout (`hfss_hammer_straps.prefix`).
* Wire bonds: `purcell_6500.set_wire_bonds(design)` + `renderer._options.update(P.WB_OPTIONS)`; short pieces:
  `ground_straps(design, 'hammer' | 'all' | 'cpw')` and `draw_hammer_straps(renderer, d, 'all')` after `render_design`.
  Bond plot: `python plot_wirebonds.py --gap <g>`.
* Fixed passes are the script defaults: `q3d/q3d_hammer_gap.py --maxp 8`, `eigenmode/hfss_hammer_straps.py --maxp 6`
  (`MESH_LINE = MESH_CAP = '10um'`), `driven/hfss_hammer_driven.py --maxp 3 --mesh-line 10um --mesh-cap 10um`.
