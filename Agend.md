# FTF Design Agenda

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
  - `my_fluxonium.py`: rectangular pads; JJ at the pocket centre; 50 µm hairpin array; U-loop flux line from the side wall; readout pad on top; charge line below, behind a ground strip.
  - `fluxonium_forprocess.py`: the same with JJ + array replaced by one lumped inductor, for simulation.
  - `my_tmon.py`: grounded coupler with flat paddles, a parallel-junction SQUID (Dolan), an L-shaped flux line, a readout pad and a charge line.
- **Numbers:** `note/design_parameters.txt` and `simulation/README.md` (latest run `tc5_x120_cl250`). Method and the full tuning log are in `note/capacitance_method.md`.
- **Simulation results:** `simulation/results/latest/` (current design) and `results/previous/` (step before); older runs only in the git history.
- **GDS settings:** `note/gds_export_settings.md`.

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

## Working agreements

- **Layout before simulation.** Show the layout first. Run simulations only after the user confirms.
- **Record every method in `note/`.** This includes the capacitance-matrix method, the κ formula and the GDS settings.
- **Q3D settings.** Medium order, 0.5 %, 6 passes for Ec / J / g. Use 12 passes at 0.2 % for κ (κ ∝ C²).
- **HFSS eigenmode.** At most 17 passes on this PC: 25 passes used more than 14 GB of RAM.
- **Fillets.** All CPW routes use 99 µm fillets. Check that every inner segment is ≥ 198 µm.
- **Resonator lengths** are measured along the filleted centre line.

## Workflow

1. **Specification:** qubit count, target Ec, resonator f and κ, readout g, control/readout port assignment.
2. **Layout:** edit the components or `ftf_layout.py`, rebuild, export, look at `fig/` (plots are made from the GDS, so fillets are real).
3. **Q3D:** use `simulation/q3d_ftf.py` (`with_resonators=True` for g). Grounded flux lines are drawn as ground in the simulation geometry, to avoid Q3D net-overlap errors.
4. **Matrix → parameters:** use `note/cap_matrix.py`.
   - Ec and J use the thesis A.14–A.15 sum/difference transform.
   - g = 2e·β·V_zpf with the λ/2 V_zpf.
   - Charge-line β and T1 come from the same matrix.
5. **κ:** run `simulation/q3d_kappa.py` (full feedline + resonators), then apply the lab formula Q_c = π / (Z0 Zr (C ω)²).
6. **Frequencies:** run `simulation/hfss_eigen_res.py`, then scale lengths by f_HFSS / f_target (`hfss_scale`).
7. **GDS:** `export_gds(..., fab=True)` gives the cheesed fabrication mask.
8. **Update** `note/design_parameters.txt`, `simulation/README.md` and this agenda; move the old `results/latest/` runs to `results/previous/`.

## Skill routing

| task | read |
|---|---|
| new / modified component, layout, routing | `skill/skill.md` |
| Q3D, HFSS, capacitance matrix, Ec / J / g / κ | `skill/skill_lom.md`, `note/capacitance_method.md` |
| physical design rules (CPW, fillets, junction placement, loss) | `skill/qubit_chip_design_guidelines.md` |
