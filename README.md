# Metal-agnet: one FTF set (fluxonium – tunable transmon coupler – fluxonium)

A Qiskit Metal layout of two fluxoniums coupled through a grounded tunable-transmon
coupler, on a 10 × 10 mm chip. It includes readout resonators, charge and flux lines,
and launch pads at fixed PCB port positions. Ec, couplings, readout g and κ are
tuned with Ansys Q3D, and the resonator frequencies are checked with HFSS eigenmode.

![chip](fig/ftf_layout_chip.png)

## Repository layout

| folder | contents |
|---|---|
| `pattern design/` | Components and the chip layout. `ftf_layout.py` (build, GDS export, plots), `my_fluxonium.py`, `fluxonium_forprocess.py` (simulation version: JJ + array → one lumped inductor), `my_tmon.py` (grounded tunable coupler), `pcbport.py` (PCB port positions) |
| `simulation/` | Ansys Q3D / HFSS scripts, sweeps in `tuning/`, results in `results/` (see [Simulation](#simulation)) |
| `note/` | `design_parameters.txt` (all final numbers), `capacitance_method.md` (method + full tuning log), `cap_matrix.py` (matrix → Ec / J / g / charge-line code), `gds_export_settings.md` (fab / cheesing GDS settings) |
| `gds/` | `ftf_layout_fab.gds` (fabrication, cheesed ground), `ftf_layout.gds` (no cheesing), `ftf_layout_forprocess.gds` (simulation geometry) |
| `fig/` | Layout plots made from the exported GDS |
| `skill/` | Design and coding rules for this workflow (Qiskit Metal components, Q3D/LOM, physical design guidelines) |
| `Agend.md` | Entry point and agenda for design sessions: status, workflow, open items |

Port use: P7 / P8 readout feedline, P1 / P6 FX flux, P2 / P5 FX charge, P3 TC flux,
P4 TC charge.

The exported GDS leaves every junction site empty. The fluxonium JJ, its 50 µm
junction array and the coupler's SQUID junctions are added separately.

## Latest design parameters (v4, CPW 10 / 6 µm)

Source: `simulation/results/latest/` — current layout, **all CPWs 10 / 6 µm (Z0 = 50.9 Ω)**.
Everything below is computed from the Q3D capacitance matrices (no retuning after the CPW
change): `cpw10_6` (FX1 + TC + FX2 with the three resonators rendered) and `kappa_cpw10_6`
(feedline + tees + resonators). The "previous" column is the 15 / 5.909 µm layout
(`simulation/results/previous/`). Full list: `note/design_parameters.txt`.

**Charging energies**, Ec = e²/2 · C̃⁻¹[i,i]

| | Ec | C_Σ | previous | target |
|---|---|---|---|---|
| FX1 | 1.310 GHz | 14.79 fF | 1.303 GHz | 1.3 GHz |
| FX2 | 1.312 GHz | 14.76 fF | 1.307 GHz | 1.3 GHz |
| TC | 0.201 GHz | 96.49 fF | 0.201 GHz | 0.2 GHz |

**Qubit–qubit coupling**, J = 4e² C̃⁻¹[i,j] / h, per charge number (n_i n_j)

| | J/2π | previous |
|---|---|---|
| FX1 – TC (f–t) | 136.9 MHz | 136.4 MHz |
| FX2 – TC (f–t) | 138.9 MHz | 138.3 MHz |
| FX1 – FX2 (f–f, direct) | 16.1 MHz | 12.1 MHz |

**Qubit–resonator coupling**, ħg = 2e β V_zpf, β = C̃⁻¹[q,r]/C̃⁻¹[r,r],
V_zpf = ω_r √(ħZ0/π) (λ/2, Z0 = 50.9 Ω, ω_r = design frequency)

| | β | V_zpf | g/2π | previous | target |
|---|---|---|---|---|---|
| FX1 – R1 (6.3 GHz) | 0.200 | 1.64 µV | 157.9 MHz (per n) | 150.9 MHz | 150 MHz |
| FX2 – R2 (6.5 GHz) | 0.193 | 1.69 µV | 157.1 MHz (per n) | 149.9 MHz | 150 MHz |
| TC – R_TC (6.0 GHz) | 0.236 | 1.56 µV | 159.4 MHz (transmon conv., EJ/Ec = 50; 178.2 MHz per n) | 149.7 MHz | 150 MHz |

g rises by √(50.9/45) ≈ 6 % with the higher Z0 (β is unchanged).

**Readout resonators** (λ/2, CoupledLineTee on one feedline P8 → P7, 99 µm fillets),
Q_c = π / (Z0 Zr (C₁ω)²), κ = f/Q_c, Z0 = Zr = 50.9 Ω, C₁ = −C[resonator, feedline]

| | f (design) | C₁ to feedline | Q_c | κ/2π | previous κ | κ target | length | tee coupling length |
|---|---|---|---|---|---|---|---|---|
| R1 (FX1) | 6.3 GHz | 11.96 fF | 5413 | 1.16 MHz | 2.08 MHz | 2 MHz | 9.268 mm | 0.689 mm |
| R_TC | 6.0 GHz | 15.50 fF | 3552 | 1.69 MHz | 3.05 MHz | 3 MHz | 9.610 mm | 0.882 mm |
| R2 (FX2) | 6.5 GHz | 11.42 fF | 5577 | 1.17 MHz | 2.08 MHz | 2 MHz | 8.983 mm | 0.662 mm |

The tee lengths were set with the 15 / 5.909 µm Cc fit (24.0 fF/mm); the 10 / 6 µm tees give
about 18.5 fF/mm, so κ is ~58 % of target. f_r is the design value: the lengths are analytic
(with the old HFSS correction) and have not been checked with HFSS for 10 / 6 µm.

**Charge (XY) lines**, Γ₁ = 8e²β²|n|²ωR/ħ, R = 50 Ω

| | C to qubit | β | T1 limit | power at chip, 10 MHz Rabi | previous T1 |
|---|---|---|---|---|---|
| FX1 | 0.322 fF | 0.0093 | 3.8 ms | −83 dBm | 2.0 ms |
| FX2 | 0.332 fF | 0.0095 | 3.6 ms | −83 dBm | 2.0 ms |
| TC | 0.038 fF | 0.00039 | 1.7 ms | −76 dBm | 0.93 ms |

FX: f01 = 0.5 GHz, |n01| = 0.1 assumed; TC: f = 5 GHz, EJ/Ec = 50 assumed.

## Running

Environment: the `qmetal` conda env (quantum-metal 0.7.4, PySide6, pyEPR, gdspy),
Ansys Electronics Desktop 2021 R2 for Q3D / HFSS (start AEDT before any simulation script).

```
cd "pattern design"
python ftf_layout.py                 # build the chip, write gds/*.gds and fig/*.png
```

In Python: `from ftf_layout import build_design; design = build_design()`. Use
`build_design(for_process=True)` for the simulation geometry, and
`export_gds(design, file, fab=True)` for the cheesed fabrication GDS.

Placing junction cells (`export_gds(..., include_junctions=True)`) needs the lab's
junction GDS library at `JJ_LIBRARY`. It is not part of this repo.

## Simulation

All scripts build the chip with `pattern design/ftf_layout.py`; the capacitance-matrix maths
is in `note/cap_matrix.py`, the method in `note/capacitance_method.md`.

| script (`simulation/`) | what it computes | output |
|---|---|---|
| `q3d_ftf.py` | Ec, J (f–t, f–f), readout g, charge-line coupling / T1 limit: Q3D of FX1 + TC + FX2 (`with_resonators=True` adds the resonators) | `<tag>_cmat_fF.csv`, `<tag>_summary.json` |
| `q3d_kappa.py` | κ: Q3D of feedline + tees + resonators, C₁ = −C[res, feedline] | `<tag>.csv`, `<tag>_cmat_fF.csv` |
| `q3d_tee.py` | tee coupling capacitance Cc vs coupling length (linear fit used by the layout) | `tee_sweep.csv` |
| `hfss_eigen_res.py` | f_r: HFSS eigenmode, resonator + tee + its qubit | `eigen_resonators.csv`, `eigen_<name>_convergence.csv` |
| `tuning/sweep_FX.py` | FX Ec / readout g around the current FX geometry | `sweep_<set>.csv` |
| `tuning/sweep_TC.py` | TC Ec / g / charge-line T1 around the current TC geometry | `sweep_<set>.csv` |
| `tuning/sweep_FX_charge_line.py` | FX charge-line tip position vs T1 limit | `sweep_fx_charge.csv` |
| `reference/analyse_fxtfx24.py` | Ec / J / g of one FTF set of FXTFX24-W7-0x-V.3.1 (comparison, `note/capacitance_method.md` §10) | – |

```
cd simulation
python -c "from q3d_ftf import Q3DRunner; r = Q3DRunner(); r.run('mytag', with_resonators=True); r.close()"
python q3d_kappa.py mytag
python q3d_tee.py 0.4 0.8 1.2          # coupling lengths in mm
python hfss_eigen_res.py FX1 TC FX2    # ~10 min per resonator
python tuning/sweep_TC.py tc5
```

New runs are written to `simulation/results/runs/`. When a run becomes the new design point,
move the current `latest/` files to `previous/` and the new ones to `latest/`, and update the
tables above and `note/design_parameters.txt`.

| folder (`simulation/results/`) | content |
|---|---|
| `latest/` | Current design (CPW 10 / 6 µm): `cpw10_6` (Ec, J, g, charge lines), `kappa_cpw10_6` (κ) |
| `previous/` | Step before (CPW 15 / 5.909 µm): `tc5_x120_cl250` (Ec, J, g, charge lines) and `sweep_tc5.csv`; `v4_pads_only` (readout-pad C used for the resonator lengths); `tee_*` (Cc fit, 15 / 5.909 µm tee); `kappa_v4d_conv` (κ); `eigen_*` (HFSS f_r) |
| `reference/` | FXTFX24 cell-33 capacitance matrix |
| `runs/` | New runs |

`*_cmat_fF.csv` is the Q3D Maxwell matrix in fF, `*_summary.json` the analysed parameters with
the geometry overrides used, `*.log` the console output, `sweep_*.csv` the table of one sweep.
Older runs are in the git history (v1–v3: commit `32bed39`; v4_final / kappa_v4c: commit `ce5674d`).

| solver | setting | use |
|---|---|---|
| Q3D | Medium order, no auto increase, 0.5 %, 6 passes (~1 min) | Ec / J / g |
| Q3D | 12 passes, 0.2 % (~2.5 min) | κ (`KAPPA_SETUP`), since κ ∝ C² |
| HFSS eigenmode | 0.05 %, ≤ 17 passes (~10 min) | f_r (25 passes used > 14 GB RAM) |
