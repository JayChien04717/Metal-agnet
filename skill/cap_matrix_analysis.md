# Skill: Analysing a capacitance matrix (Q3D → Ec, J, g, β, κ)

How to turn a Q3D capacitance matrix of a Qiskit Metal layout into qubit / coupler / resonator
parameters, and how to check that the matrix can be trusted. This file is the procedure.
The numbers of the current chip (targets, results, tuning log) are in
`note/capacitance_method.md` and `note/design_parameters.txt`; the code is `note/cap_matrix.py`
and `simulation/q3d_ftf.py`.

---

## 1. Workflow

1. **Simulation geometry.** Use the `for_process=True` version of each component: junction + array
   replaced by one lumped inductor, same pads / pocket / lines (`fluxonium_forprocess.py`,
   `my_tmon.py` with `flux_as_ground=True`).
2. **Choose what is rendered.** Qubits + coupler for Ec / J. Add the resonators for g (the readout
   node must be the whole resonator, not a small floating pad). Feedline + tees + resonators
   (no qubits) for κ.
3. **Open-terminate every pin** so readout pads and charge lines become their own conductors.
4. **Run Q3D** (settings in section 6) → Maxwell matrix `C` in fF, saved as `<tag>_cmat_fF.csv`.
5. **Pick the nodes** (section 2), transform, invert, read off parameters (sections 3–5).
6. **Sanity-check** the matrix (section 7) *before* quoting any number.
7. **Record** in `note/design_parameters.txt` / `note/capacitance_method.md`.

---

## 2. Which conductors are nodes

Q3D net names are `<polygon>_<component>` (`pad_top_FX1`, `tmon_TC`, `readout_line_FX1`,
`charge_line_FX1`, `ground_main_plane`).

| conductor | treat as |
|---|---|
| qubit pads / islands | node |
| readout resonator (pad + rendered line, one net) | node (only when you want g or the resonator load) |
| charge line | node only for β / T1 of the drive line; otherwise grounded |
| flux lines (shorted to ground) | ground; draw them as ground plane in the sim geometry |
| ground plane | reference, always dropped |

**Rule:** a conductor that is not in the node list is grounded, i.e. its row and column are
*deleted* from the Maxwell matrix. This is correct for lines terminated in 50 Ω at the qubit
frequency, and wrong for floating metal – floating metal must be a node.

---

## 3. Maxwell matrix and sign conventions

* Q3D returns the **Maxwell** matrix: `C_ii > 0` is the *total* capacitance of node i to everything
  (all other nodes + ground), `C_ij < 0` for i ≠ j.
* Mutual (coupling) capacitance between i and j is `Cc_ij = −C_ij`.
* Capacitance of node i to ground = row sum `Σ_j C_ij` (ground is the reference).
* Ec of a single isolated island is `e²/(2 C_ii)` only if all other nodes are grounded; with several
  nodes use the inverse matrix (next section).

---

## 4. Ec and J (multi-node)

For a fluxonium / transmon with two pads, use the sum / difference coordinate so that only the
differential (qubit) mode carries charge number `n`:

```python
C   = cmat.loc[nodes, nodes].to_numpy() * 1e-15            # F
M   = block_diag([[1, 1], [1, -1]] per pad pair, 1 elsewhere)
Mi  = np.linalg.inv(M)
Ct  = Mi.T @ C @ Mi                                        # C~
Ci  = np.linalg.inv(Ct)                                    # C~^-1
Ec_i = e**2 / 2 * Ci[i, i] / h                             # i = difference coordinate of the qubit
J_ij = 4 * e**2 * Ci[i, j] / h                             # coupling per n_i n_j
```

* The inverse is essential: `Ci[i,i]` already contains the loading of every other node.
* Target total capacitance for a given Ec: `C_Σ = e² / (2 h Ec)` (Ec = 1.3 GHz → 14.9 fF,
  Ec = 0.2 GHz → 96.9 fF).
* A **ground-referenced transmon / coupler island** is one coordinate (row of M = 1).
* J depends on the node set. Floating readout pads that are not connected to the real resonator add
  a fake coupling path, so quote J from the run **with resonators**.

---

## 5. Readout g, β, charge line, κ

**Voltage-division ratio** (insensitive to the node's own self capacitance):

```python
beta = Ci[q, r] / Ci[r, r]                 # = Cc / C_Σ,q
```

**Readout coupling** (λ/2 resonator, open ends):

```python
V_zpf = w_r * np.sqrt(hbar * Z0 / np.pi)
g     = 2 * e * beta * V_zpf / h           # per charge number n (matrix element excluded for fluxonium)
g_tr  = g * (32 * Ec / EJ) ** 0.25         # transmon convention
```

* `w_r` comes from the design (or HFSS eigenmode), *not* from `1/(2 Z0 C_line)` of a coarse Q3D
  mesh – that overestimates f_r by ~30 %.
* Multiply g by the qubit charge matrix element `|<0|n|1>|` for the transition you read out.
* To lower g: increase the pad–pad distance (`readout_line_options.pad_sep`) or shrink the pad.

**Charge line (drive port)**, same β with the line as node:

```python
gamma1 = 8 * e**2 * beta**2 * n01**2 * w * R / hbar    # T1 limit from the port, R = 50 Ohm, T = 0
V_rabi = h * f_rabi / (2 * e * abs(beta) * n01)         # amplitude for a given Rabi rate
```

Check T1 limit ≫ intrinsic T1 target (design for ≥ ms for the drive port) and that the required
drive amplitude is reachable.

**Resonator linewidth κ.** Q3D of feedline + tees + resonators (qubits left out):

```python
C_1 = -C[res, feedline]                         # coupling capacitance, F
Q_c = np.pi / (Z0 * Zr * (C_1 * w) ** 2)
kappa = f / Q_c                                  # κ/2π
```

Inverse for design: `C_c = sqrt(π / (Q_c Z0 Zr)) / ω`. κ ∝ C_1², so a ±0.3 fF mesh error is ±3 %.
Use more passes (section 6).

---

## 6. Q3D settings that work

```python
lom.sim.setup.freq_ghz = 5.0
lom.sim.setup.max_passes = 6           # 12 for κ
lom.sim.setup.min_passes = 2
lom.sim.setup.percent_error = 0.5      # 0.2 for κ
lom.sim.setup.solution_order = "Medium"
lom.sim.setup.auto_increase_solution_order = False
```

* ~1 min for the chip. `High` order + auto increase + 0.1 % + 20 passes took > 15 min and looked hung.
* Noise at these settings: ≈ 1 % in Ec / J / g.
* **Convergence check:** run once with more passes / tighter error; Ec must move < 1–2 %.
* Use a fresh Q3D design name per run (stale NULL bodies after a failed render).

---

## 7. Sanity checks before trusting the matrix

| check | expected | if it fails |
|---|---|---|
| matrix symmetric | `C_ij = C_ji` within solver noise (< 1 %) | mesh not converged |
| signs | diagonal > 0, off-diagonal < 0 | wrong net grouping |
| row sum ≥ 0 | = capacitance to ground | missing / merged net |
| pad self capacitance | sum of the pad's row ≈ expected from geometry (parallel-plate / CPW estimate within ~2×) | wrong unit (fF vs F) or geometry unit |
| `C_Σ` vs target | Ec within the tolerance of the design target | geometry or node set |
| pass-to-pass | last two passes agree < 1 % | add passes |
| grounded vs floating | nothing that should float is missing from the node list | add node |
| resonator loading | Q3D `C_pad` used in the length formula, ≈ 1 % in f_r per halving of the pad | – |

Independent cross-check: estimate the pad–pad / pad–ground capacitance with a simple model
(coplanar strip / CPW formula using εeff = (1+εr)/2) and compare; a large mismatch (factor ≳ 2)
points to a units or net problem, not physics. Also compare with the previous design's matrix
(`results/previous/`): a small geometry change should give a small matrix change.

---

## 8. Renderer / Ansys pitfalls

* Launch AEDT (`ansysedt.exe`) before any script. COM check:
  `win32com.client.Dispatch('Ansoft.ElectronicsDesktop').GetAppDesktop().GetVersion()`.
  Re-register stale COM with `reg_ansysedt.exe`.
* Stack: `pandas==2.2.3` (pyEPR breaks on pandas 3); `delim_whitespace=True` → `sep=r"\s+"`.
* Avoid long, tiny, U-shaped or heavily filleted `path` / `path_cut` geometry in Q3D
  (`Failed to sweep the object`, `start and end point coincident`, `Subtract` failures). Convert to
  `poly`, reduce fillets, or leave it out of a capacitance-only model.
* Zero-area slivers from `draw.subtract` → drop them in the component (`_clean`), else
  *Cannot create rectangle with zero length/width*.
* `Net … overlap`: a flux line shorted to ground becomes a second net touching ground → draw it as
  ground.
* `capacitance_all_passes` may hold numpy arrays; wrap in `pd.DataFrame(m, index=…, columns=…)`.
* To find the failing body: monkeypatch `QAnsysRenderer.render_element_path` to print component /
  name / subtract / length; or subtract the cutouts one at a time.

---

## 9. Reporting checklist

- [ ] Which node set and which components were rendered (with / without resonators)
- [ ] Q3D settings and number of passes / final ΔC
- [ ] Ec, J, g, β (and T1 limit) with the formula used
- [ ] Matrix saved (`results/runs/<tag>_cmat_fF.csv`) and sanity checks of section 7 passed
- [ ] Number quoted from the run *with resonators* where relevant
