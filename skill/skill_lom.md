# LOM / Q3D Capacitance Matrix Workflow

Use this note when extracting Qiskit Metal capacitance matrices with Ansys Q3D and converting them into LOM qubit/coupler/resonator parameters.

## Ansys / COM

- Qiskit Metal's Ansys renderer uses pyEPR / pywin32 COM automation. In this desktop environment, Ansys COM may fail inside the sandbox with `server execution failed`; rerun Ansys automation commands with elevated execution.
- Minimal COM check:

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

- Reference `fluxonium_v2_withJJ.py` for Q3D-friendly fluxonium geometry:
  - Main capacitor pads and arms are `poly`.
  - Charge/readout/flux connector shapes are also mostly `poly`.
  - Junctions and superinductor elements are placed in the `junction` table.
  - Material is usually not manually assigned in the component. Q3D renderer assigns metal as thin conductor/PEC and chip material comes from `design.chips.main.material`, e.g. `"silicon"`.
- Avoid sending long, tiny, U-shaped, or heavily filleted `path` / `path_cut` geometries into Q3D. AEDT ACIS can fail with:
  - `Specified start and end point coincident`
  - `Failed to sweep the object`
  - `Cannot calculate normal vector, divide by zero`
  - `Subtract` failures on ground cutouts.
- For FTF, the stable Q3D model kept capacitive pads / readout / charge / coupler and excluded:
  - Flux-bias pickup paths.
  - External PCB feed routes.
  - Fluxonium array/superinductor path and path-cut geometries.
- Use a fresh Q3D design name for each run to avoid stale NULL bodies after failed AEDT renders.

## Debugging Renderer Failures

- If `SweepAlongPath` fails, monkeypatch `QAnsysRenderer.render_element_path` to print component, qgeometry name, subtract flag, and length before calling the original method.
- If `subtract_from_ground()` fails, temporarily subtract one shape at a time and print the failing body name.
- Once the failing geometry is known, either convert it to `poly` geometry following `fluxonium_v2_withJJ.py`, reduce/remove fillets, or exclude it from the capacitance-only Q3D model.

## Matrix Extraction

- Use `LOManalysis(design, "q3d")`.
- For capacitance extraction only:

```python
lom.sim.setup.name = "FTF_LOM_Q3D"
lom.sim.setup.freq_ghz = 5.0
lom.sim.setup.max_passes = 15
lom.sim.setup.min_passes = 2
lom.sim.setup.min_converged_passes = 2
lom.sim.setup.percent_error = 0.5
lom.sim.setup.percent_refinement = 30

lom.sim.run(
    name=f"FTF_LOM_clean_{timestamp}",
    components=selected_components,
    open_terminations=open_pins,
    box_plus_buffer=False,
)
cap_matrix = lom.sim.capacitance_matrix
```

- `capacitance_all_passes` may contain numpy arrays instead of DataFrames; convert with `pd.DataFrame(pass_matrix, index=cap_matrix.index, columns=cap_matrix.columns)` before saving.

## LOM Transform

For the FTF calculation, use the 8-node order:

```python
ORDER = [
    "pad_left_FL",
    "pad_right_FL",
    "cross_TC",
    "pad_left_FR",
    "pad_right_FR",
    "readout_launch_FL",
    "prime_cpw_LR_tee",
    "readout_launch_FR",
]
```

and transformed labels:

```python
NEW_LABELS = [
    "FL_sum",
    "FL_diff",
    "TC",
    "FR_sum",
    "FR_diff",
    "res_FL",
    "res_c",
    "res_FR",
]
```

Use:

```python
C = cap_ff.loc[ORDER, ORDER].to_numpy(dtype=float) * 1e-15
C_tilde = np.linalg.inv(M) @ C @ np.linalg.inv(M)
C_tilde_inv = np.linalg.inv(C_tilde)
```

with:

```python
M = np.array([
    [1,  1, 0, 0,  0, 0, 0, 0],
    [1, -1, 0, 0,  0, 0, 0, 0],
    [0,  0, 1, 0,  0, 0, 0, 0],
    [0,  0, 0, 1,  1, 0, 0, 0],
    [0,  0, 0, 1, -1, 0, 0, 0],
    [0,  0, 0, 0,  0, 1, 0, 0],
    [0,  0, 0, 0,  0, 0, 1, 0],
    [0,  0, 0, 0,  0, 0, 0, 1],
], dtype=float)
```

## Parameter Formulas

Using zero-based Python indices on `C_tilde_inv`:

```python
EC1 = 0.5 * e**2 * C_tilde_inv[1, 1]
ECc = 0.5 * e**2 * C_tilde_inv[2, 2]
EC2 = 0.5 * e**2 * C_tilde_inv[4, 4]

J1c = 4 * e**2 * C_tilde_inv[1, 2]
J2c = 4 * e**2 * C_tilde_inv[2, 4]
J12 = 4 * e**2 * C_tilde_inv[1, 4]
```

For fluxonium-resonator coupling in the convention where the resonator matrix element is included but the qubit matrix element is excluded:

```python
hbar_g = e * np.sqrt(np.pi * hbar / (2 * Z0)) * C_tilde_inv[q_index, r_index]
g_over_2pi_hz = hbar_g / h
```

Use `Z0 = 50 ohm` for the lambda/4 resonator special case unless another resonator impedance is known.

For transmon-style convention, include the qubit charge matrix element:

```python
hbar_g = e * np.sqrt(np.pi * hbar / (2 * Z0)) * (32 * EC / EJ)**0.25 * C_tilde_inv[q_index, r_index]
```

## Charge-Line Couplings

For voltage-source / charge-line couplings, keep the charge/control line nodes as the `y` block and the qubit/resonator/coupler nodes as the `x` block. With the same variable transform:

```python
C_tilde = M_inv @ Cxx @ M_inv
C_tilde_inv = np.linalg.inv(C_tilde)
Cxy_tilde = Cxy @ M_inv
beta = Cxy_tilde @ C_tilde_inv
```

The Hamiltonian drive term is:

```python
H_drive = 2 * e * V_line * beta[line, coordinate] * n_coordinate
```

The dimensionless `beta` is often the most useful charge-line coupling number. For a voltage scale, report:

```python
drive_per_mV_over_h_GHz = 2 * e * beta * 1e-3 / h / 1e9
```

For the current FTF matrix, charge-line rows are:

```python
CHARGE_LINES = [
    "charge_launch_FL",
    "charge_line_TC",
    "charge_launch_FR",
]
```

and the relevant coordinates are usually `FL_diff`, `TC`, and `FR_diff`.

## Current FTF Output Files

- `FTF_capacitance_matrix.csv`: raw Q3D Maxwell capacitance matrix in fF.
- `FTF_capacitance_matrix_table.csv`: publication-style table with an Index row.
- `FTF_Ctilde_fF.csv`: transformed capacitance matrix.
- `FTF_Ctilde_inverse_1_per_F.csv`: inverse transformed capacitance matrix.
- `FTF_qubit_parameters_from_LOM.csv`: charging energies and coupling strengths.
- `scripts/calc_ftf_lom_parameters.py`: reproducible calculation script.
- `scripts/run_ftf_q3d_lom_debug.py`: notebook-driven Ansys Q3D debug runner.
