# LOM / Q3D Capacitance Matrix Workflow

Use this note when extracting Qiskit Metal capacitance matrices with Ansys Q3D and converting them into LOM qubit/coupler/resonator parameters.

## Ansys / COM

- Qiskit Metal's Ansys renderer uses pyEPR / pywin32 COM automation. In this desktop environment, Ansys COM may fail inside the sandbox with `server execution failed`; rerun Ansys automation commands with elevated execution.
- Launch AEDT (`ansysedt.exe`) before any script. Minimal COM check:

```powershell
& C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe -c "import win32com.client; app=win32com.client.Dispatch('Ansoft.ElectronicsDesktop'); d=app.GetAppDesktop(); print(d.GetVersion())"
```

- If COM registration is stale, re-register AEDT 2021.2:

```powershell
& 'C:\Program Files\AnsysEM\AnsysEM21.2\Win64\reg_ansysedt.exe'
```

## Package Compatibility Fixes

- Old pyEPR uses pandas APIs that break on pandas 3. Prefer `pandas==2.2.3` for this stack.
- If pyEPR calls `pd.read_csv(..., delim_whitespace=True)`, replace with `sep=r"\s+"` or keep pandas below 3.
- Qiskit Metal LOM may assume pandas Series for `chi_in_MHz` / `gbus`; support numpy arrays with:

```python
lambda x: x.iloc[0] if hasattr(x, "iloc") else x[0]
```

- pyEPR variation indexing may need `int(variation)` instead of `ureg(variation)` when indexing tuple/list variations.

## Geometry Lessons

- Simulation geometry: `pattern design/fluxonium_forprocess.py` (JJ + array → one lumped
  inductor, same pads / pocket / lines) and `my_tmon.py` with `flux_as_ground=True`.
  `build_design(for_process=True)` builds the chip from them.
- Capacitor pads, readout pads, charge-line tips and the coupler are `poly`; junctions and the
  lumped inductor are in the `junction` table. Material comes from
  `design.chips.main.material` (`silicon`); metal is rendered as a thin conductor.
- Avoid sending long, tiny, U-shaped, or heavily filleted `path` / `path_cut` geometries into
  Q3D. AEDT ACIS can fail with:
  - `Specified start and end point coincident`
  - `Failed to sweep the object`
  - `Cannot calculate normal vector, divide by zero`
  - `Subtract` failures on ground cutouts.
- Grounded flux lines are drawn as part of the ground plane in the simulation geometry;
  otherwise Q3D reports net-overlap errors.
- Boolean ops can leave zero-area slivers; drop them in the component (`_clean` in
  `my_tmon.py`), otherwise `CreateRectangle` fails.
- Use a fresh Q3D design name for each run to avoid stale NULL bodies after failed AEDT renders.

## Debugging Renderer Failures

- If `SweepAlongPath` fails, monkeypatch `QAnsysRenderer.render_element_path` to print component, qgeometry name, subtract flag, and length before calling the original method.
- If `subtract_from_ground()` fails, temporarily subtract one shape at a time and print the failing body name.
- Once the failing geometry is known, either convert it to `poly` geometry (as in `fluxonium_forprocess.py`), reduce/remove fillets, or exclude it from the capacitance-only Q3D model.

## Matrix Extraction

`simulation/q3d_ftf.py` wraps this; the settings it uses:

```python
lom = LOManalysis(design, "q3d")
lom.sim.setup.freq_ghz = 5.0
lom.sim.setup.max_passes = 6            # 12 for kappa (KAPPA_SETUP)
lom.sim.setup.min_passes = 2
lom.sim.setup.percent_error = 0.5       # 0.2 for kappa
lom.sim.setup.solution_order = "Medium"
lom.sim.setup.auto_increase_solution_order = False

lom.sim.run(components=components, open_terminations=open_pins)
cap_matrix = lom.sim.capacitance_matrix   # Maxwell matrix, saved as <tag>_cmat_fF.csv
```

- `capacitance_all_passes` may contain numpy arrays instead of DataFrames; convert with
  `pd.DataFrame(pass_matrix, index=cap_matrix.index, columns=cap_matrix.columns)` before saving.

## Node names (Q3D nets = `<polygon>_<component>`)

| node | meaning |
|---|---|
| `pad_bot_FX1`, `pad_top_FX1` | FX1 pads → sum / difference coordinate |
| `tmon_TC` | coupler island |
| `pad_top_FX2`, `pad_bot_FX2` | FX2 pads |
| `readout_line_FX1`, `readout_line_FX2`, `readout_line_TC` | resonator nets (readout pad + rendered resonator) |
| `charge_line_FX1`, `charge_line_FX2`, `charge_line_TC` | charge lines |
| `ground_main_plane` | ground |

Every conductor not in the node list is treated as grounded (its rows / columns are dropped).

## Parameter Formulas (`note/cap_matrix.py`)

Sum / difference transform of each FX pad pair (thesis A.14):

```python
C = cmat.loc[nodes, nodes].to_numpy() * 1e-15
C_tilde = Mi.T @ C @ Mi                # Mi = inv(M), M rows (1, 1) / (1, -1) per pad pair
Ci = np.linalg.inv(C_tilde)
Ec_i = e**2 / 2 * Ci[i, i]             # FX difference coordinate, coupler
J_ij = 4 * e**2 * Ci[i, j] / h         # per charge number n_i n_j
```

Readout g (`readout_g`), with the resonator rendered (λ/2, distributed):

```python
beta  = Ci[q, r] / Ci[r, r]                    # voltage division, = Cc / C_sigma
V_zpf = w_r * np.sqrt(hbar * Z0 / np.pi)       # lambda/2, Z0 = 45 Ohm, w_r = design frequency
g     = 2 * e * beta * V_zpf / h               # per n, fluxonium matrix element excluded
g_tr  = g * (32 * Ec / EJ) ** 0.25             # transmon convention (EJ/Ec = 50 assumed)
```

Do not use the lumped formula `e sqrt(pi hbar / 2 Z0) C~^-1[q, r]` with the whole resonator
net as one node (`g_FX1-R1` in the summary, ≈ 67 MHz): it does not apply to a λ/2 resonator.

Charge line (`charge_line`), same β with the line as node:

```python
gamma1 = 8 * e**2 * beta**2 * n01**2 * w * R / hbar   # T1 limit, R = 50 Ohm, T = 0
V_rabi = h * f_rabi / (2 * e * abs(beta) * n01)        # drive amplitude for a Rabi rate
```

Linewidth (`simulation/q3d_kappa.py`): C_1 = −C[resonator, feedline] from Q3D of feedline +
tees + resonators, Q_c = π / (Z0 Zr (C_1 ω)²), κ = f / Q_c, Z0 = Zr = 45 Ω.

## FTF scripts in this repo

* `simulation/q3d_ftf.py` — `Q3DRunner().run(tag, with_resonators=True)`: renders FX1, TC, FX2 (+ resonators), saves `results/runs/<tag>_cmat_fF.csv` / `_summary.json`, prints Ec, J, g and charge-line numbers via `note/cap_matrix.py`.
* `simulation/q3d_kappa.py` — Q3D of feedline + tees + resonators; C_1 = −C[res, feedline]; κ from the lab formula.
* `simulation/q3d_tee.py` — one CoupledLineTee: Cc vs coupling length (linear fit used in the layout).
* `simulation/hfss_eigen_res.py` — HFSS eigenmode per resonator (resonator + tee + its qubit).
* `simulation/tuning/` — sweeps around the current geometry (`sweep_TC.py`, `sweep_FX.py`, `sweep_FX_charge_line.py`).
* Results: `simulation/results/latest/` (current design), `results/previous/` (step before), new runs in `results/runs/`; see `simulation/README.md`.
* Full method and tuning log: `note/capacitance_method.md`; final table: `note/design_parameters.txt`.

## Lessons from the FTF v4 runs

* **Launch AEDT first** (`ansysedt.exe`); the renderer connects over COM. Licenses are logged in `%TEMP%\.ansys\ansyscl.*.log` — useful to see whether a "stuck" run is really solving.
* **Solver settings:** Q3D with `solution_order='High'`, auto increase, 0.1 % and 20 passes took > 15 min and looked hung. Medium order, no auto increase, 0.5 %, 6 passes runs in ~1 min and is enough for Ec / J / g (≈1 % noise). κ ∝ C_1², so use 12 passes / 0.2 % for κ (still ±0.3 fF ≈ ±3 % in κ).
* **Net overlap errors** (`Net 'flux_line_FX1' and Net 'flux_bias_line_TC' overlap`): flux lines shorted to ground become separate nets touching the ground. Draw them as ground in the simulation geometry and do not open-terminate them.
* **Zero-area polygons** from boolean ops break `CreateRectangle`; clean them in the component.
* **Readout g needs the resonator in the model.** With only the small floating readout pad as node, C~^-1[q,r] (and A.18) is far too large. Render the resonator (it merges with the pad into one net) or use g = 2e·β·V_zpf with β = C~^-1[q,r]/C~^-1[r,r] (insensitive to the node's own C) and V_zpf = ω_r·sqrt(ħ Z0/π) for λ/2. Take ω_r from the design, not from 1/(2 Z0 C_line) of a coarse Q3D mesh.
* **J depends on which nodes are kept:** small floating readout pads add a spurious coupling path; quote J from the run with resonators.
* **Resonator loading:** use the Q3D self capacitance of the readout pad (`C_pad`) in the length formula; halving the pad changed it from 39 to 24 fF (≈ 1 % in f_r).
* **HFSS eigenmode:** keep qubits far below the resonators (fluxonium lumped 200 nH, coupler `Lj` 50 nH) so the resonator is only dispersively pulled; ≤ 17 passes (25 passes → hf3d > 14 GB RAM); expect ±10 MHz.
