# FTF Design Agenda (Design3: readout pads for g ~ 80 MHz, see README)

Last update: 2026-10-21

Read this file first at the start of a design session. Then read the `skill/` file
that matches the task (see *Skill routing*).

## Current state (v4)

One FTF set (FX1 – TC – FX2) on a 10 × 10 mm chip. Everything lives in this repo:

- **Layout:** `pattern design/ftf_layout.py`.
  - All tuned dimensions are in `FX_OPTIONS`, `FX_PER_QUBIT`, `TMON_OPTIONS` and `RESONATORS`.
  - `build_design()` builds the chip.
  - `export_gds()` writes `gds/`.
- **Components:**
  - `my_fluxonium.py`: rectangular pads; JJ at the pocket centre (2 µm gap); 55 × 10 µm hairpin array (2 µm wires and leads); U-loop flux line from the side wall, tip 15 µm from the array end; readout pad on top; charge line below, behind a ground strip.
  - `fluxonium.py` (`FluxoniumArrayJJ`): the simulation model for **every** run (Q3D via `build_design(for_process=True)`, HFSS eigenmode, mutual inductance): loop = lead line from the stem ends with the JJ as its short middle gap (qubit centre); U = two lumped-L sides (100 nH each = 200 nH) joined by a metal closing end; JJ = port in the mutual run.
  - `my_tmon.py`: grounded coupler with flat paddles, a parallel-junction SQUID (Dolan), an L-shaped flux line, a readout pad and a charge line.
- **Numbers:** `design_parameters.txt` and `README.md` (latest runs `fxarray_fl15`, `kappa_cpw10_6`). Method and the full tuning log are in `../Design1/notes/capacitance_method.md`.
- **Simulation results:** `simulation/results/latest/` (current design) and `results/previous/` (step before); older runs only in an earlier repository (commits `32bed39` / `ce5674d`, not available here).
- **GDS:** junction-free, `../skill/gds_export.md`.

Targets (all met):

| | target |
|---|---|
| Ec_FX | 1.3 GHz |
| Ec_TC | 0.2 GHz |
| readout g FX / TC | 150 MHz |
| resonators | 6.3 / 6.0 / 6.5 GHz, κ 2 / 3 / 2 MHz |

## Open items

- [x] CPW changed to 10 / 6 um (50.9 Ohm): Q3D re-simulated (`cpw10_6`, `kappa_cpw10_6`). Ec unchanged; g ~158 MHz (+5 %); kappa 1.16 / 1.69 / 1.17 MHz (~58 % of target).
- [ ] Retune for 10 / 6 um: tee coupling lengths for kappa (Cc ~18.5 fF/mm), readout gaps for g = 150 MHz, then HFSS f_r.
- [x] Coupler charge line on the right of the stem (x = +120 um), L-shaped flux line shifted left: Q3D re-simulated, T1 limit 0.93 ms, Ec / g unchanged.
- [ ] HFSS eigenmode on the corrected resonator lengths (9.286 / 9.647 / 9.000 mm). Expected on target within ±10 MHz.
- [ ] Junction / array cells: the GDS leaves the junction sites empty. Add the real JJ, the 50 µm array and the SQUID cells (cell name `jj_array_50um` for the array).
- [ ] Charge-line T1 / drive estimates use an assumed f01 = 0.5 GHz and |n01| = 0.1. Redo them with the real EJ, EL.
- [x] Flux-line mutual inductance (HFSS driven modal, `simulation/driven/hfss_mutual_fx.py`): M = 0.44 pH at fl_sep 20 um, I_Phi0 = 4.7 mA (`Design1/notes/mutual_inductance_method.md`).
- [x] 2026-11-03: real JJ array (55 x 10 um, 2 um) in the layout, flux tip 15 um: M = 0.575 pH (I_Phi0 = 3.6 mA); Q3D `fxarray_fl15`: Ec_FX 1.322 / 1.326 GHz, g_FX 161 / 160 MHz (all within ~2 % of `cpw10_6`).
- [x] 2026-11-03: readout pad -> 10 x 40 um stub, gap 10 / 11.45 um (FX1 / FX2): Q3D `ro_w10` g 80.9 / 80.5 MHz (was 161 / 160), Ec_FX 1.340 / 1.342 GHz, chi ~ 1.8 MHz (`Design1/notes/zz_method.md`). Layout updated; resonator lengths / f_r (C_line -1 %) and kappa still to retune.
- [ ] Check that the `junction_1` cell fits the new 2 um JJ gap (`jj_length` in FX_OPTIONS).

## Working agreements

- **Layout before simulation.** Show the layout first. Run simulations only after the user confirms.
- **Record every method in `note/`.** This includes the capacitance-matrix method, the κ formula and the GDS settings.
- **Q3D settings.** Medium order, 0.5 %, 6 passes for Ec / J / g. Use 12 passes at 0.2 % for κ (κ ∝ C²).
- **HFSS eigenmode.** At most 17 passes on this PC: 25 passes used more than 14 GB of RAM.
- **Fillets.** All CPW routes use 99 µm fillets. Check that every inner segment is ≥ 198 µm.
- **Resonator lengths** are measured along the filleted centre line.

## Workflow

1. **Specification:** qubit count, target Ec, resonator f and κ, readout g, control/readout port assignment.
2. **Layout:** edit the components or `ftf_layout.py`, rebuild, export, look at `fig/layout/` (plots are made from the GDS, so fillets are real).
3. **Q3D:** use `simulation/q3d/q3d_ftf.py` (`with_resonators=True` for g). Grounded flux lines are drawn as ground in the simulation geometry, to avoid Q3D net-overlap errors.
4. **Matrix → parameters:** use `note/cap_matrix.py`.
   - Ec and J use the thesis A.14–A.15 sum/difference transform.
   - g = 2e·β·V_zpf with the λ/2 V_zpf.
   - Charge-line β and T1 come from the same matrix.
5. **κ:** run `simulation/q3d/q3d_kappa.py` (full feedline + resonators), then apply the lab formula Q_c = π / (Z0 Zr (C ω)²).
6. **Frequencies:** run `simulation/eigenmode/hfss_eigen_res.py`, then scale lengths by f_HFSS / f_target (`hfss_scale`).
7. **GDS:** `export_gds(build_design(jj_metal=False), file, fab=True)` gives the junction-free cheesed mask (`../skill/gds_export.md`).
8. **Update** `design_parameters.txt`, `README.md` and this agenda; move the old `results/latest/` runs to `results/previous/`.

## Skill routing

| task | read |
|---|---|
| new / modified component, layout, routing | `skill/skill.md` |
| Q3D, HFSS, capacitance matrix, Ec / J / g / κ | `skill/cap_matrix_analysis.md`, `Design1/notes/capacitance_method.md`; flux-line M: `skill/mutual_inductance_analysis.md`, `Design1/notes/mutual_inductance_method.md` |
| layout DRC / physical design rules (CPW, fillets, junction placement, loss) | `skill/chip_drc.md` |
