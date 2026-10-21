# simulation/

Ansys Q3D / HFSS runs for the FTF design (fluxonium FX1 – transmon coupler TC – fluxonium FX2),
driven through qiskit-metal. Start Ansys Electronics Desktop 2021 R2 before running any script.
All scripts build the chip with `pattern design/ftf_layout.py`; the capacitance-matrix maths is
in `note/cap_matrix.py`, the method in `note/capacitance_method.md`, the full parameter list in
`note/design_parameters.txt`.

## Latest design parameters

Source: `results/latest/` (current layout: TC charge line at x = +120 µm, tip 250 µm below the bar,
TC flux line shifted 25 µm left). Q3D runs with the three resonators rendered.

**Charging energies** (`tc5_x120_cl250`), Ec = e²/2 · C̃⁻¹[i,i]

| | Ec | C_Σ | target |
|---|---|---|---|
| FX1 | 1.303 GHz | 14.86 fF | 1.3 GHz |
| FX2 | 1.307 GHz | 14.83 fF | 1.3 GHz |
| TC | 0.201 GHz | 96.44 fF | 0.2 GHz |

**Qubit–qubit coupling** (`tc5_x120_cl250`), J = 4e² C̃⁻¹[i,j] / h, per charge number (n_i n_j)

| | J/2π |
|---|---|
| FX1 – TC (f–t) | 136.4 MHz |
| FX2 – TC (f–t) | 138.3 MHz |
| FX1 – FX2 (f–f, direct) | 12.1 MHz |

**Qubit–resonator coupling** (`tc5_x120_cl250`), ħg = 2e β V_zpf, β = C̃⁻¹[q,r]/C̃⁻¹[r,r],
V_zpf = ω_r √(ħZ0/π) (λ/2, Z0 = 45 Ω, ω_r = design frequency)

| | β | g/2π | target |
|---|---|---|---|
| FX1 – R1 (6.3 GHz) | 0.203 | 150.9 MHz (per n) | 150 MHz |
| FX2 – R2 (6.5 GHz) | 0.195 | 149.9 MHz (per n) | 150 MHz |
| TC – R_TC (6.0 GHz) | 0.236 | 149.7 MHz (transmon conv., EJ/Ec = 50; 167.3 MHz per n) | 150 MHz |

Use `readout_FX1/FX2/TC → g/2pi` in the summary json. The `g_FX1-R1` key (≈ 67 MHz) is the
lumped-J formula with the whole resonator as one node and does not apply to a λ/2 resonator.

**Readout resonators** (λ/2, CoupledLineTee on one feedline P8 → P7, 99 µm fillets)

| | f target | f HFSS¹ | κ/2π² | κ target | C₁ to feedline | Q_c | length | tee coupling length |
|---|---|---|---|---|---|---|---|---|
| R1 (FX1) | 6.300 GHz | 6.289 GHz | 2.08 MHz | 2 MHz | 18.10 fF | 3023 | 9.286 mm | 0.775 mm |
| R_TC | 6.000 GHz | 6.026 GHz | 3.05 MHz | 3 MHz | 23.56 fF | 1966 | 9.647 mm | 0.995 mm |
| R2 (FX2) | 6.500 GHz | 6.489 GHz | 2.08 MHz | 2 MHz | 17.23 fF | 3132 | 9.000 mm | 0.743 mm |

¹ `eigen_resonators.csv`, simulated at 9.302 / 9.606 / 9.015 mm; lengths then scaled by
f_HFSS/f_target (not re-simulated), so expected f_r ≈ 6.30 / 6.00 / 6.50 GHz (±10 MHz).
² `kappa_v4d_conv.csv`: Q_c = π / (Z0 Zr (C₁ω)²), κ = f/Q_c, Z0 = Zr = 45 Ω (±3 %).

**Charge (XY) lines** (`tc5_x120_cl250`), Γ₁ = 8e²β²|n|²ωR/ħ, R = 50 Ω

| | C to qubit | β | T1 limit | power at chip, 10 MHz Rabi |
|---|---|---|---|---|
| FX1 / FX2 | 0.447 fF | 0.0128 | 2.0 ms (f01 = 0.5 GHz, \|n01\| = 0.1) | −86 dBm |
| TC | 0.051 fF | 0.00053 | 0.93 ms (f = 5 GHz, EJ/Ec = 50) | −79 dBm |

## Scripts

| script | what it computes | output |
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
python q3d_ftf.py                      # pads only, tag 'baseline'
python -c "from q3d_ftf import Q3DRunner; r = Q3DRunner(); r.run('mytag', with_resonators=True); r.close()"
python q3d_kappa.py mytag
python q3d_tee.py 0.4 0.8 1.2          # coupling lengths in mm
python hfss_eigen_res.py FX1 TC FX2    # ~10 min per resonator
python tuning/sweep_TC.py tc5
```

New runs are written to `results/runs/`. When a run becomes the new design point, move the
current `latest/` files to `previous/` and the new ones to `latest/`.

## results/

| folder | content |
|---|---|
| `results/latest/` | Current design. `tc5_x120_cl250` (Ec, J, g, charge lines) and its sweep `sweep_tc5.csv` / `tc5.log`; `v4_pads_only` (readout-pad C used for the resonator lengths); `tee_*` (Cc fit); `kappa_v4d_conv` (κ); `eigen_*` (HFSS f_r) |
| `results/previous/` | Step before: `v4_final` (TC charge line still left of the stem; FX numbers the same within 0.5 %), `kappa_v4c_conv` (κ run that set `C_extra` in the layout) |
| `results/reference/` | FXTFX24 cell-33 capacitance matrix |
| `results/runs/` | New runs |

`*_cmat_fF.csv` is the Q3D Maxwell matrix in fF, `*_summary.json` the analysed parameters with
the geometry overrides used, `*.log` the console output, `sweep_*.csv` the table of one sweep.
Older runs (v1–v3 tuning) are in the git history (commit `32bed39` on GitHub).

## Solver settings

| solver | setting | use |
|---|---|---|
| Q3D | Medium order, no auto increase, 0.5 %, 6 passes (~1 min) | Ec / J / g |
| Q3D | 12 passes, 0.2 % (~2.5 min) | κ (`KAPPA_SETUP`), since κ ∝ C² |
| HFSS eigenmode | 0.05 %, ≤ 17 passes (~10 min) | f_r (25 passes used > 14 GB RAM) |
