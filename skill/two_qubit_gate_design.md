# Skill: Two-qubit gate design (transmon – tunable coupler – transmon, flux CZ)

Operating procedure for selecting parameters and verifying a flux-activated CZ gate on two flux-tunable transmons (Q1, Q2) coupled through a flux-tunable transmon coupler (TC).

* **Physics, Hamiltonian, and formulas:** [note/tunable_coupler_cz_method.md](../note/tunable_coupler_cz_method.md).
* **Design-specific scripts, values, and runs:** each Design's `simulation/qubit/`, `notes/`, and `design.json`. (Move reusable CZ modules to `code/qubit/` once shared by multiple Designs.)
* **Environment:** `scqenv` (`scqubits`) for quantum simulations; EM inputs come from `qmetal` runs ([em_simulation.md](em_simulation.md)).

---

## 1. Workflow

1. **Extract Ec and g from EM:** use a core Q3D run (islands, coupler, claws, XY lines; flux lines drawn as ground, [cap_matrix_analysis.md](cap_matrix_analysis.md)) or an HFSS eigenmode fit.
   - Record the upstream EM run and pass count (pass counts are not convergence evidence).
   - State accepted EM uncertainty (≈ 3 % on Ec and g is typically sufficient for gate-design decisions unless stricter accuracy is requested).
2. **Take junctions from the layout:** use the Design's layout `LJ` inputs.
   - EJ is a user design input—**never retarget it automatically** ([parameter_provenance.md](parameter_provenance.md)).
   - Label any re-solved junction set (`retarget`) explicitly as a *proposal*; never write proposals back into the layout or `design.json`.
3. **Static scan:** with qubits at their sweet spots, sweep coupler flux and record dressed f1, f2, fc, α1, α2, and ZZ vs. coupler frequency.
4. **Choose idle (ZZ-off) point and CZ scheme** (Section 3). Record the adopted scheme and user decision source.
5. **Run time-domain dynamics:** scan gate time T (coarse grid, then refine around optimum) and record conditional phase, leakage, coherent fidelity F, flat-top frequencies, whether amplitude hits its bound, and peak populations of |20⟩, |02⟩, and coupler-excited states.
6. **Record via `code/general/run_record.py`:** set `upstream` to the EM run, log all `parameters` (levels, `dt`, bounds, `dissipation`), and tag `status_of_inputs` (`target` / `proposal` / `result`).
7. **Update Design records:** add a row to `results/index.csv`, summarize in the Design's CZ note, and mark earlier-geometry results `historical`.

---

## 2. Units and convention pitfalls (check before comparing numbers)

| Quantity | Convention | Pitfall |
|---|---|---|
| `g_ij` | **Charge units** (`g_ij n_i n_j`) | Not exchange coupling `g_ex σ⁺σ⁻` (which is ~1.5–2× larger; see [method note §1](../note/tunable_coupler_cz_method.md#1-hamiltonian)). Never compare without converting. |
| `ZZ` | `E11 − E10 − E01 + E00` (dressed) | Sign depends on regime; always read sign and zero from the full diagonalised model. |
| Conditional phase `φ` | `(φ11 − φ10 − φ01 + φ00) mod 360°` | With `ZZ > 0`, `φ` counts *down* from 360°—a value of 200–300° can mean "too little ZZ", not an optimiser failure. Verify with a direct amplitude scan. |
| Fidelity `F` | Computational-subspace average gate fidelity, virtual Z removed | Coherent only unless T1/T2 are included—always state explicitly. |

---

## 3. Choosing the scheme

* **Coupler-only CZ** (qubits fixed at sweet spots; only TC pulsed): requires a real ZZ-off point—i.e., **straddling regime** (`|Δ12| < |α|`), `g12` cancelling coupler-mediated coupling at idle (verify `g12` with Q3D), and `f_c,max` well above the qubits without colliding with readout frequencies. Gate time is bounded by `T_flat ≳ 1 / (2 |ZZ_on|)` at the coupler proximity bound (default `f_c ≥ f_Q + 0.3 GHz`).
* **Q1 + TC CZ** (Q1 pulsed to `f_Q1 = f_Q2 − α1` for `|11⟩–|20⟩` plus TC pulse): works for `|Δ12| > |α|` and gives shorter gates. Pulse Q1 *towards* resonance without crossing Q2 (avoids `|01⟩–|10⟩` swap). Not usable when `f_Q2 − α1` lies above Q1's sweet spot.

---

## 4. Checklist before quoting a gate

- [ ] Idle ZZ value, sign, and bias method (exact ZZ zero vs. minimum `|ZZ|`).
- [ ] Coupler idle frequency checked against readout/filter frequencies and harmonics (≥ 100 MHz clearance).
- [ ] Optimal T is smooth across neighbouring T values (not an isolated point or optimiser miss).
- [ ] Stated whether amplitude hits the coupler bound (if yes, time limit is `ZZ_on`, not pulse shape).
- [ ] Peak populations checked: `|11⟩–|20⟩` hybridisation is expected; coupler-excited population > ~1 % means the bound is too close.
- [ ] Leakage quoted at chosen T and at ±10 % (timing-sensitive near `|20⟩`).
- [ ] Explicitly stated: coherent only (no T1/T2), truncation levels & `dt`, whether pulse shape is unoptimised, and whether idle flux-noise sensitivity was checked.
- [ ] Upstream EM revision matches the layout revision; otherwise mark `needs_revalidation`.

---

## 5. Not covered yet

Pulse-shape optimisation (Slepian / Martinis–Geller), open-system simulation (T1, T2, flux noise), and XY-drive crosstalk to the coupler during single-qubit gates. Extend Sections 1 and 4 and the method note when added.
