# Design4: transmission Purcell filter, 6.7 GHz (interdigital finger couplers)

Last update: 2026-11-04. Python: see `skill/skill.md` section 9 (conda env `qmetal`).

## Current state (layout only, NOT simulated)
- `pattern design/interdigital_cap.py`: `InterdigitalCap` (N fingers, finger width / gap, overlap, end gap, bus, CPW pins a / b).
- `pattern design/purcell_layout.py`: P1 -> feed -> cap -> lambda/2 meander -> cap -> feed -> P6, 10 x 10 mm chip, CPW 10 / 6 um, 99 um fillet.
- Spec assumed: f0 = 6.7 GHz, kappa_F = 200 MHz (Q_F 33.5). Change `F0`, `KAPPA_F` at the top of the script.
- Numbers (analytic): Cc 101 fF per port -> N = 10 fingers, 6 / 6 um, overlap 199.1 um (Alley formula); lambda/2 length 7.486 mm (7.406 mm filleted between cap pins).
- GDS `gds/purcell_filter.gds`, figures `fig/layout/`.

## Open items
- [ ] User confirms layout / bandwidth (kappa_F) / finger count.
- [ ] Q3D of one interdigital cap -> real Cc (replace Alley), retune `OVERLAP`.
- [ ] HFSS driven (S21) or eigenmode: f0 and kappa_F; set `HFSS_SCALE`.
- [ ] Add readout resonators coupled to the filter (not in this device yet).
- [ ] Faint horizontal hairlines in the GDS plot at y ~ 1.2 / 1.7 / 2.0 mm: check in KLayout whether they are real slivers.

## Working agreements
Same as Design3: layout before simulation, methods in `../note/`, 99 um fillets, lengths on the filleted centre line.

## Results
- 2026-11-04 Q3D `cap_v1` (cap_in alone, taper 100 um, 10 fingers 6 / 6 um, overlap 199.1 um; 12 passes): Cc = 97.99 fF (target 101.1, Alley 101.0),
  Q_c 71.3, kappa_F = 188 MHz (target 200). C to ground per side 77.4 fF (large: low-impedance comb region, matching is the open issue).
  Next: overlap ~ 205 um for 101 fF (Cc roughly linear in overlap, not re-run), then HFSS driven S21 / S11.
  Files: simulation/results/cap_v1_cmat_fF.csv, simulation/cap_v1.log.

## Status 2026-11-05
- Kept as is (user decision): series transmission filter (`purcell_layout.py`, standalone, s21_v3: f0 6.73 GHz, 3 dB BW 210 MHz, S11 -22 dB)
  and the Design3 + filter version (`design3_purcell.py`: P8 -> cap -> filter line with 3 tees -> cap -> P7, resonators 6.6 / 6.7 / 6.8 GHz).
- Combined-version HFSS sweeps disagree (Interpolating comb_v1: peak 6.647 GHz, no dips; Fast comb_v2: broad peak ~6.55 GHz, one sharp feature at 6.74 GHz);
  Discrete check `comb_v3_discrete` (61 points, 10 MHz) decides which one to trust. Not resolved yet.
- Reference topology (Heinsoo 2018, arXiv:1801.07904): per qubit lambda/4 readout resonator -(C_c at x_c)- individual lambda/4 Purcell filter -(finger cap)- common feedline
  with input cap C_in = 40 fF; kappa_P 40 MHz, J 10 MHz, kappa_R = 1/2 (kappa_P - Re sqrt(-16 J^2 + (kappa_P - 2 i Delta)^2)), 160 MHz spacing. Not built.
