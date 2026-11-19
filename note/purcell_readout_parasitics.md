# Purcell filter + readout resonator: layout for hammer-only coupling

General method for the shared module [`PurcellReadout`](../pattern%20design/resonator/purcell_readout.py).
The goal is that **frequency comes from length** and **filter–resonator coupling comes from the hammer C_j only**.
Then J follows from one small Q3D coupon, and f_p, f_r follow from the filleted length plus Q3D end loads.
Design-specific numbers and runs are in [Design8](../Design8/README.md).

## 1. Where unwanted coupling comes from

Earlier Designs showed three paths besides the hammer (Design4 notes, [method note §5](purcell_filter_method.md)):

| Path | Seen in | Size |
|---|---|---|
| Side-by-side CPW runs (capacitive + inductive) | Design4 at 200 µm, ~1.5 mm of parallel run | 6.0 MHz of 12.0 MHz HFSS J; the hammer gave 2.9 MHz |
| Ground cut by stub + hammer pocket (return current) | Design4 at 250 µm without straps | J 68 MHz → 23 MHz with 2 straps → 20 MHz with straps on every CPW |
| Meander self/neighbour coupling at < 400 µm | skill.md §11 | 1–2 MHz |

## 2. Layout rules in the module

1. **Opposite-side routing.** Each main line passes straight through its T. The filter's open arm goes to C_p on
   one side; its shorted arm turns away from the resonator (jog, then a meander on the far side). The resonator is
   the point mirror. The only side-by-side run is the T region: 2 × (`tee_half` + `lead`) = 300 µm at the default
   settings.
2. **Line spacing ≥ 400 µm** at the T (skill.md §11). The ground strip between the lines is 378 µm wide (10/6 µm
   CPW) and is cut only by the stubs and the hammer pocket.
3. **Meanders ≥ 400 µm from everything.** `jog` (400 µm) puts the first meander run 800 µm from the other line.
   The run pitch is 400 µm.
4. **Ground straps / air bridges at the hammer** (`strap_points()`): one across each stub, halfway between the T
   cut and the pocket, and one across each main line on both sides of each T (Wallraff C_AB). That is six bridges
   per module. In HFSS they are PEC straps; on the mask they are air bridges. Every other CPW gets qiskit auto bonds
   every 350 µm (skill.md §13).
5. **Perpendicular pieces only near the hammer.** The stubs run perpendicular to the main lines. Jogs and meander
   steps that face the other line are perpendicular or end-on (exempt in §11).

Alternatives that were considered:
- **Perpendicular main lines** (filter T rotated 90°): no parallel run at all, but the hammer bar must bend, so
  the coupon is no longer a straight two-pin piece.
- **Wider spacing:** at 600 µm, c′_m halves again (table below), but the longer stubs move the hammer away from
  the lumped-element regime.
- **Ground fingers between the T regions:** not needed. The 378 µm strip is already solid ground. What matters is
  that it stays one piece, which the straps ensure.

## 3. Numbers

**Distributed coupling (2-D extraction).** These are cross-section values from
[`code/resonator/cpw_pair_cap2d.py`](../code/resonator/cpw_pair_cap2d.py)
([data](purcell_readout_data/cpw_pair_cap2d.json)). The solver is a finite-volume Laplace solve with
zero-thickness metal, a 500 µm substrate and natural outer boundaries. As a check, it gives 0.546 fF/mm for a
16/9 µm pair at 200 µm (ε_r 11.5), against 0.563 fF/mm from the Design4 Q3D pair coupons. Its c′_m changes by less
than 1.3 % between grid steps of 1 µm and 0.25 µm. The single-line self C′ converges slowly, about 5 % high at
the finest grid, so take C′ from the conformal value.

| Centre spacing (µm) | 150 | 200 | 300 | 400 | 600 | 800 |
|---|---|---|---|---|---|---|
| c′_m 10/6 µm, ε_r 11.9 (fF/mm) | 0.400 | 0.229 | 0.106 | 0.062 | 0.028 | 0.015 |
| c′_m 16/9 µm, ε_r 11.9 (fF/mm) | 0.990 | 0.564 | 0.260 | 0.151 | 0.068 | 0.036 |

k_L ≈ k_C (for example 3.3e-4 against 3.6e-4 at 400 µm, 10/6 µm). The inductive part is therefore as large as the
capacitive part and cannot be neglected.

**Design8 module (10/6 µm, 400 µm, f 6.5 GHz).** Q3D run
[`20261118_232405_q3d_module_b3d273d0`](../Design8/results/resonator/20261118_232405_q3d_module_b3d273d0/summary.json)
used 8 passes. No case reported solver convergence; last-pass changes were 0.8–4.9 %.

| Quantity | Value | Source |
|---|---|---|
| Hammer C_j (coupon mutual − straight-line reference mutual) | 1.971 fF | Q3D (fit to Design4 data: 1.908 fF) |
| Full-module static mutual (ends rendered open) | 2.076 fF | Q3D |
| Static non-hammer mutual | 0.105 fF (5.3 % of C_j) | Q3D difference |
| Same, side-by-side runs only (2-D c′_m × run lengths) | 0.037 fF | 2-D model |
| J from the hammer (closed form, Q3D C_j) | 13.19 MHz | model |
| Parallel-run J bound (\|J_C\| + \|J_L\|) | ≤ 0.30 MHz (2.2 %) | model |

The static non-hammer mutual (0.105 fF) is about 3× the 2-D side-by-side estimate. The rest is near-field coupling
around the T regions and pocket, plus Q3D mesh noise: the straight reference alone has 0.032 fF of mutual, and its
two lines differ by 2.6 fF to ground although they are symmetric. Static mutual capacitance also overstates J,
because most of the extra part sits where the mode voltages are small.

**Design8 HFSS check** (open fixture, lengths from the Q3D-loaded model, nothing fitted; [Design8](../Design8/README.md#2-validation-results-open-fixture-revision-2026-11-18_r1)):

| | Model | HFSS driven, 6 passes | HFSS eigenmode, 6 passes |
|---|---|---|---|
| J | 13.19 MHz (hammer only) | 13.39 MHz (×1.02) | 16.71 MHz (×1.27) |
| f_p / f_r | 6.500 / 6.500 GHz | +8.2 / −9.4 MHz | +10.0 / +15.8 MHz |
| κ_p | 26.8 MHz | 25.6 MHz | 47.2 MHz |

The driven J matching the hammer-only value within 2 % is the direct evidence that the layout removes the other
coupling paths. In Design4 (old pair) the same comparison gave HFSS 12.0 MHz against 2.9 MHz from the hammer.
The bound of 0.30 MHz from the side-by-side runs (2.2 %) is consistent with the residual.

DESIGN56_PLACEHOLDER

## 4. What is not simulated

- The qubit-loaded channel. The lengths for the fluxonium and transmon chips use the Design3 Q3D pad capacitances,
  which were extracted for the old 15/5.9 µm CPW. Run a qubit-loaded Q3D before using those lengths.
- Air-bridge loading on the mask. HFSS uses thin PEC wires, not real bridges.
- The perpendicular-routing variant, and line spacings other than 400 µm, in Q3D or HFSS.
- Q3D convergence beyond 8 passes.
