# Skill: Two-qubit gate design (transmon - tunable coupler - transmon, flux CZ)

This file is the operating procedure for choosing parameters and checking a flux-activated CZ on two flux-tunable transmons Q1, Q2 coupled through a flux-tunable transmon coupler TC.
- **Physics and formulas:** [`note/tunable_coupler_cz_method.md`](../note/tunable_coupler_cz_method.md).
- **Numbers for a chip:** that Design's `notes/` and `design.json`. Never put a chip's values in this file.

The current implementation is in Design6. It should move to `code/qubit/` when a second Design uses it.

| Entry point | Does |
|---|---|
| `Design6/simulation/qubit/cz_scqubits.py` | `Model` (Hamiltonian), `load_params`, `evolve`, `analyse_U`; legacy CLI for the Q1 + TC scheme. Its `main` writes fixed paths, so call its functions from a run-recorded wrapper. |
| `Design6/simulation/qubit/cz_coupler_flux.py` | coupler-only scheme: `idle_point`, `cz_at` |
| `Design6/simulation/qubit/cz_r9.py` | run-recorded wrapper. It builds the parameters from a Q3D core run plus the layout `LJ` and runs both schemes. Copy this pattern. |

Environment: `scqenv` (scqubits). Q3D and HFSS inputs come from `qmetal` runs.

---

## 1. Workflow

1. **Get Ec and g from EM.** Use a core Q3D (islands, coupler, claws, XY lines; flux lines drawn as ground), following [`skill/cap_matrix_analysis.md`](cap_matrix_analysis.md), or an HFSS eigenmode fit.
   - Record the source run and the pass. Pass counts are not convergence evidence.
   - State the accepted error: about 3 % on Ec and g is enough for gate-design decisions unless the user asks for more.
2. **Get the junctions from the layout.** Use the fab inputs (`LJ` in the Design layout).
   - EJ is a user design input. **Do not retarget it automatically.**
   - A parameter set with re-solved junctions (`retarget`) is a *proposal*. Label it so, and never write it back to the layout or `design.json`.
3. **Static scan.** With the qubits at their sweet spots, sweep the coupler flux. Record:
   - the dressed f1, f2, fc and α1, α2;
   - ZZ vs coupler frequency.
4. **Choose the idle (ZZ-off) point and the scheme** (section 3).
5. **Run the dynamics.** For each candidate:
   - scan the gate time T (a coarse grid, then refine near the best point);
   - record the conditional phase, leakage, F, the coupler and qubit frequencies at the flat top, whether the amplitude is at its bound, and the peak populations of |20⟩, |02⟩ and the coupler-excited states.
6. **Make a new run** with `code/general/run_record.py`.
   - `upstream` is the EM run.
   - `parameters` holds the sets, levels, dt, bounds and `dissipation`.
   - `status_of_inputs` says whether the inputs are target, proposal or result.
7. **Write up.** Use the Design's CZ note and its parameter overview, and add a row to `results/index.csv`. Mark earlier-geometry results `historical`.

## 2. Units and conventions (check before comparing any number)

| Quantity | Convention used here | Pitfall |
|---|---|---|
| g_ij | **charge units**: H term g n_i n_j, with g = 4e² (C⁻¹)_ij / h | This is not the exchange coupling g σ⁺σ⁻. g_ex = g_n ⟨0\|n_i\|1⟩⟨0\|n_j\|1⟩, with ⟨0\|n\|1⟩ ≈ (EJ/8Ec)^(1/4)/√2 ≈ 1.2–1.5 for EJ/Ec ≈ 50–150. The exchange value is therefore about 1.5–2× *larger*. Never compare across conventions without converting. |
| ZZ | E11 − E10 − E01 + E00 of the dressed levels | Its sign depends on the regime (section 3). |
| Conditional phase | φ = (φ11 − φ10 − φ01 + φ00) mod 360° | With ZZ > 0 the phase counts *down* from 360°. A result of 200–300° can therefore mean "too little ZZ". It does not mean the optimiser failed. Check with a direct amplitude scan. |
| F | average gate fidelity, F = (\|Tr(U_CZ† U)\|² + Tr(U†U)) / 20 on the computational subspace, virtual Z removed | Coherent only unless T1 and T2 are added. Always say so. |
| Asymmetric SQUID | d = (EJa − EJb)/(EJa + EJb) | The flux range and the sweet spot follow from d. |

## 3. Choosing the scheme

**Coupler-only CZ** (qubits fixed at their sweet spots, only the coupler pulsed) needs a real ZZ-off point:
- **Straddling regime:** |Δ12| < |α|. When |Δ12| > |α|, ZZ keeps one sign at every coupler bias and never crosses zero. The best you get is a minimum |ZZ|.
- **g12 matched:** the direct coupling g12 must cancel the coupler-mediated coupling at the idle bias. The off-point position depends steeply on g12, so it needs its own Q3D check.
- **High coupler f_max:** the coupler must sit far enough above the qubits at idle. Check it against the readout frequencies.

**Gate-time limit:**
- The conditional phase is ∫ZZ dt, so T_flat ≳ 1/(2 |ZZ_on|), with ZZ_on taken at the coupler bound.
- The bound is f_c ≥ f_Q + margin; 0.3 GHz is the default, and closer gives more coupler excitation.
- Longer or shorter ramps do not beat this limit.

**Q1 + TC CZ** (Q1 pulsed to the |11⟩–|20⟩ condition, f_Q1 = f_Q2 − α1, plus a coupler pulse):
- It works for |Δ12| > |α| and gives the shortest gates.
- Pulse Q1 *towards* the resonance without crossing Q2. Crossing Q2 causes a |01⟩–|10⟩ swap error.
- It does not apply when f_Q2 − α1 lies above Q1's sweet spot (straddling sets), because Q1 cannot move up.

**Decision record.** Write down which scheme is adopted and why. The adoption is the user's decision; record its source.

## 4. Checks before quoting a gate

- [ ] Idle ZZ and its sign.
- [ ] Idle bias method: ZZ zero or minimum |ZZ|.
- [ ] Coupler frequency at idle vs readout and filter frequencies, including harmonics or collisions within about 100 MHz.
- [ ] The best gate time is not an isolated point: neighbouring T values behave smoothly, or the failure is shown to be an optimiser miss.
- [ ] Amplitude at the bound? If yes, the time limit is ZZ_on, not the pulse.
- [ ] Peak populations: |11⟩–|20⟩ hybridisation is expected. Coupler-excited population above about 1 % means the bound is too close.
- [ ] Leakage is timing-sensitive near |20⟩. Quote the leakage at the chosen T and at ±10 %.
- [ ] State: coherent only (no T1/T2); levels (e.g. 4 × 5 × 4) and dt; that the pulse shape is not optimised (if true); and that flux-noise sensitivity at idle is not checked (if true).
- [ ] The EM source revision matches the layout revision of the gate result. Otherwise mark it `needs_revalidation`.

## 5. Not covered yet

- Pulse-shape optimisation (Slepian / Martinis–Geller).
- Open-system simulation (T1, T2, flux noise).
- XY-drive crosstalk to the coupler during single-qubit gates.

When any of these is added, extend sections 1 and 4 and the method note.
