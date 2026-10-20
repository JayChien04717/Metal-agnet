# simulation/

Ansys Q3D / HFSS runs for the FTF design, driven through qiskit-metal.
Start Ansys Electronics Desktop 2021 R2 before running any script. All of them
build the chip with `pattern design/ftf_layout.py`. The capacitance-matrix maths
is in `note/cap_matrix.py`; the method is in `note/capacitance_method.md`.

## Main scripts (the ones that produce the final numbers)

| script | what it computes | how | output |
|---|---|---|---|
| `q3d_ftf.py` | Ec (FX, TC), couplings J (f–t, f–f), readout g (f–r, t–r), charge-line coupling / T1 limit | Q3D of FX1 + TC + FX2 (`with_resonators=True` adds the three resonators), then `cap_matrix` | `<tag>_cmat_fF.csv`, `<tag>_summary.json` |
| `q3d_kappa.py` | resonator linewidth κ | Q3D of feedline + tees + resonators; C_1 = −C[res, feedline] → Q_c = π / (Z0 Zr (C_1 ω)²) | `<tag>.csv`, `<tag>_cmat_fF.csv` |
| `q3d_tee.py` | tee coupling capacitance Cc vs coupling length | Q3D of one CoupledLineTee | `tee_sweep.csv` (linear fit used by the layout) |
| `hfss_eigen_res.py` | resonator frequencies f_r | HFSS eigenmode: resonator + tee + its qubit, one design per resonator | `eigen_resonators.csv`, `eigen_<name>_convergence.csv` |

Usage:

```
python q3d_ftf.py                      # pads only, tag 'baseline'
python -c "from q3d_ftf import Q3DRunner; r = Q3DRunner(); r.run('mytag', with_resonators=True); r.close()"
python q3d_kappa.py mytag
python q3d_tee.py 0.4 0.8 1.2          # coupling lengths in mm
python hfss_eigen_res.py FX1 TC FX2    # ~10 min per resonator
```

New runs are written to **`results/runs/`**.

## tuning/: sweeps used to tune the design

Each script runs a named set of geometry variants and writes a table
`sweep_<set>.csv` to `results/runs/`.

| script | tuned | knobs |
|---|---|---|
| `tuning/sweep_v1_Ec.py` | FX Ec → 1.3 GHz, TC Ec → 0.2 GHz (first design) | FX pad size, pocket, FX spacing, coupler stem / bar |
| `tuning/sweep_TC.py` | TC Ec, TC readout g = 150 MHz, TC charge-line T1 | stem length, readout pad width, charge-line tip |
| `tuning/sweep_FX.py` | FX Ec, FX readout g = 150 MHz | FX pad size, readout pad gap / height |
| `tuning/sweep_FX_charge_line.py` | FX charge-line coupling vs T1 limit | charge-line tip position `cl_sep` |

## reference/

| script | what |
|---|---|
| `reference/analyse_fxtfx24.py` | Ec / J / g of one FTF set of the FXTFX24-W7-0x-V.3.1 design, from that notebook's saved Q3D matrix (comparison in `note/capacitance_method.md` §10) |

## results/

| folder | content |
|---|---|
| `results/final/` | **The runs behind every number in `note/design_parameters.txt`.** `v4_final` (Ec, J, g, charge lines), `v4_pads_only` (readout-pad capacitance used for resonator lengths), `tee_*` / `tee_sweep.csv` (Cc fit), `kappa_v4c_conv` (C_extra), `kappa_v4d_conv` (final κ), `eigen_*` (HFSS f_r), `tc5_x120_cl250` (final TC charge line / flux line geometry: TC Ec, g, charge-line T1) |
| `results/history/1_v1_Ec_tuning/` | First design: pad / pocket / FX spacing sweeps (`baseline`, `s1_*`–`s3_*`, `final_default`) |
| `results/history/2_v1_readout_g/` | First g with resonators rendered (`final_with_res`) |
| `results/history/3_v2_hairpin_array/` | JJ at the centre + 50 µm hairpin array, side U flux line (`v2_with_res`) |
| `results/history/4_v3_resonators_TC_readout/` | Three resonators + TC readout; TC stem / pad / charge-line tuning (`v3_full`, `tc1`–`tc3`, `v3_final`) |
| `results/history/5_v3_FX_charge_line/` | FX charge-line strength sweep (`fxcl_*`, `v3b_final`) |
| `results/history/6_v4_pocket_readout/` | Smaller FX pocket, half-width readout pad, g = 150 MHz tuning (`v4_base`, `fx1_*`, `v4_tuned`) |
| `results/history/7_kappa_6pass/` | κ runs with the 6-pass mesh, too noisy for κ (`kappa_v4`, `v4b`, `v4c`) |
| `results/history/8_hfss_first_15pass/` | First HFSS run on FX1 (15 passes, not converged) |
| `results/history/9_v4_TC_lines/` | TC charge line moved right + tip sweep (`tc4_*` at x = 160 µm, `tc5_*` at x = 120 µm) |
| `results/reference/` | FXTFX24 cell-33 capacitance matrix |
| `results/runs/` | New runs |

File types: `*_cmat_fF.csv` is the Q3D Maxwell matrix in fF, `*_summary.json` the
analysed parameters with the geometry overrides used, `*.log` the console output,
and `sweep_*.csv` the table of one sweep.

## Settings

| solver | setting | use |
|---|---|---|
| Q3D | Medium order, no auto increase, 0.5 %, 6 passes (~1 min) | Ec / J / g |
| Q3D | 12 passes, 0.2 % (~2.5 min) | κ (`KAPPA_SETUP`), since κ ∝ C² |
| HFSS eigenmode | 0.05 %, ≤ 17 passes (~10 min) | f_r. 25 passes used more than 14 GB of RAM. |
