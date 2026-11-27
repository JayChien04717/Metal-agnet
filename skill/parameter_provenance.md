# Parameter provenance, run storage, and confirmation rules

Applies to new Designs, new runs, and organization of older results. Do not overwrite historical raw data; add provenance through sidecars without inventing missing evidence.

## 1. Roles and states

- **role:** `target`, `input`, `derived` (calculated/fitted), or `result`. Downstream inputs must retain upstream provenance.
- **status:** `draft`, `recorded` (stored does not mean validated), `validated` (supported by evidence), `unknown`, or `not_applicable`.
- **applicability:** `current`, `historical`, `needs_revalidation`, or `unknown`. Track validation status separately from geometry applicability.

## 2. Source boundaries

| Parameter | Source and requirements |
|---|---|
| EJ, EL | User settings or quantum-model solutions under explicit targets/constraints; record whether fixed or adjustable |
| Lumped Lj / Larray | EM inputs; preserve the conversion direction to/from EJ/EL |
| Capacitance matrix | Raw Q3D result; record net order, reference ground, units, and elimination |
| EC, charge coupling | Matrix derivation or EM circuit fitting; do not label as direct EM outputs |
| Linear modes, S/Z | HFSS; record boundaries, ports, mesh, frequency range, and convergence |
| Qubit f01, α, ZZ | Quantum model; record EJ/EL/EC, bias, and truncation |
| Fidelity, leakage | Dynamics calculation; record pulse, time step, algorithm, dissipation, and definitions |

- HFSS linear frequencies are not nonlinear qubit f01. If `scqubits` solves for EJ using HFSS-fitted EC, record it as `"scqubits target solve, dependent on HFSS-fitted EC"`, not as EJ obtained from HFSS.
- Never automatically retarget user-fixed EJ/EL. Mark unclear historical authorization `unknown`. If a transmon model has no independent EL, record `not_applicable` rather than `0`.

## 3. Per-parameter records

Required fields: `value`, `unit`, `role`, `source_kind`, `source_ref` (exact file and field/function), `method`, `depends_on`, `constraints`, `status`, `applicability`, `geometry_revision`, and `run_id`. Unknown values may be `null` with an explanation.

Explicitly label E/h in GHz. Distinguish f from ω and κ from κ/2π. Specify the Hamiltonian coefficient definition for coupling (charge vs. exchange). Φ/Φ0 is dimensionless.

### 3.1 Standard three-tier parameter presentation (README and delivery PDF)

Both `Design<N>/README.md` and per-design sheets in fabrication report PDFs (`skill/fabrication_output.md` §5) must separate parameters into **three explicit tiers** so fabrication engineers and circuit designers can read their exact numbers without cross-referencing mixed tables:

1. **Tier 1 — Fabrication Team Targets (`fab_junction_targets`):**
   - Only what is needed to fabricate and inspect the chip: chip/die ID, revision, source fab GDS + SHA256, marker / test-junction / logo summary, and per-device junction targets.
   - **Transmon / SQUID:** `E_J1`, `E_J2` (GHz, `E/h`), symmetry (`symmetric: EJ1 = EJ2` or asymmetric ratio / `d`), `L_J1 / L_J2` (nH), `R_n1 / R_n2` (kΩ, state Δ), and exact junction-site centres `(x, y)` + orientation + SQUID spacing when the e-beam site draws the junctions.
   - **Fluxonium:** small-junction `E_J` and superinductance `E_L` (GHz, `E/h`), realisation (array JJ count `N`, per-JJ `E_J,arr`), `R_n` targets, and mask widths when drawn.
   - **Test junctions:** count, type, coordinates, and target `R_n` / `E_J`.
   - Missing fab values must print **`NOT SET`**; never substitute HFSS lumped inductances.

2. **Tier 2 — EM & Circuit Design Parameters (separated by subsystem):**
   - **2A. Per-qubit / coupler:** for each element (`Q1`, `Q2`, `TC`, …), `E_C` (MHz), `C_Σ` (fF), idle bare & dressed `f_01` (GHz), anharmonicity `α` (MHz), and operating flux `Φ/Φ_0`.
   - **2B. Qubit–qubit & qubit–coupler (`Q–Q`, `Q–TC`) couplings:** charge-basis `g_n` (MHz, `g_n n_i n_j`), exchange `g_ex` (MHz, `g_ex σ⁺σ⁻`), mutual capacitance `C_ij` (fF), and idle `ZZ` (kHz).
   - **2C. Qubit–readout (`Q–R`) couplings:** per qubit `C_qr` (fF), `g_qr` (MHz), dispersive shift `2χ` (MHz), `|2χ|/κ_r,eff`, and critical photon number `n_crit`.
   - **2D. Purcell–readout (`P–R`) & resonator modes:** per channel `f_r` (GHz), `f_p` (GHz), detuning `f_r − f_p` (MHz), `J` (MHz), `κ_p` (`κ_in / κ_out`, MHz), effective linewidth `κ_r,eff` (MHz), and Purcell `T_1` (state model: Sete vs. Lorentzian).
   - **2E. Flux-line & control (when applicable):** per line `|M|` (pH), `I_Φ0` (mA), gate current `I_gate` (mA), `L_fl` / SQUID `L_geo` (pH), and crosstalk matrix `M_ij / M_jj` (%).

3. **Tier 3 — Time-Domain `scqubits` & Coherence Analysis (when available):**
   - **Two-qubit gate dynamics:** gate scheme, pulse shape & ramp time, optimal `T_gate` (ns), conditional phase (deg), coherent fidelity `F_coh`, leakage, and any `T_1`-limited estimate (explicitly stating assumed `T_1` and that it is an analytic estimate, not a master-equation simulation).
   - **Readout & dephasing estimates:** single-shot `F_ro` & integration time `τ` (with assumed `T_1`, `n̄`, `η`), and Clerk–Utami photon shot-noise `T_φ` vs. `n_th` / `T_eff` (using `κ_r,eff`).
   - If time-domain calculations have not been run for the current geometry, explicitly state `Not run` or `needs_revalidation` rather than silently omitting or substituting older data.

## 4. Directory structure and per-run records

- **Shared tools:** `code/general/`, `code/qubit/`, `code/resonator/`, `code/flux/`.
- **Per-Design categories:** `layout/`, `simulation/{general,qubit,resonator,flux}/`, and matching `results/{general,qubit,resonator,flux}/`. Do not add `jobs/`, `support/`, `release/`, or `sweep/` layers.
- **New runs (`code/general/run_record.py`):** under a unique run ID, store `run.json`, `code/`, `raw/`, and `summaries/` / `figures/` / `logs/` (plus `inputs/` when input snapshots are needed). Record time/timezone, command, working directory, environment/versions, geometry revision, commit/dirty state, input hashes, upstream runs, execution status, and validation evidence.
- **Scoped, flat `code/` snapshots:** `<run>/code/` must snapshot **only** the entrypoint and the project Python modules it actually imports (via AST + runtime `sys.modules`), stored flat under `code/` (parent-prefixed only on basename collisions), with original repo-relative `path`, `snapshot`, and `sha256` recorded in `run.json`. Never mirror the entire repository tree or copy unrelated Design folders into `<run>/code/` (avoids bloat and Windows `MAX_PATH` failures). Pass non-Python inputs (e.g., `geometry.json`, GDS) explicitly via `inputs=`.
- A commit hash alone does not identify a dirty workspace: save actual input snapshots.
- **Re-analysis:** create a separate run referencing `upstream` instead of overwriting adopted data.
- **Legacy data:** keep old mixed results in `general/legacy/` without rewriting original records (older migrated runs may still contain `manifest.json` / `analysis/`). Connect legacy entry points that write fixed filenames to a new run output directory before execution; a migration index is not historical provenance evidence.
- **Solver side files:** never write `hfss_eig_f_convergence.csv` or similar files to the repository root (`project_paths.configure()` and `run_record.create_run()` redirect them to `raw/`).

## 5. Adopting results

- Select `adopted_results` explicitly in `design.json` rather than picking the newest filename. README presents rounded summaries linked to exact source values.
- Mark earlier results `needs_revalidation` after changes to geometry, junction models, or boundaries. Restore `current` only with a rerun or concrete applicability check.

## 6. Final version confirmation before commit/push

- Before committing or pushing, prepare a reviewable summary of the final candidate, parameters, adopted runs, GDS/figures, and unvalidated items; obtain the user's explicit confirmation of the final version.
- Reuse explicit confirmation while version and parameters remain unchanged; any adjustment invalidates it and requires reconfirmation.
- Version confirmation does not itself authorize `git commit` / `git push`; perform Git operations only when requested by the user.
- README parameters and results must match only the final confirmed version. Keep older values in `notes/` or historical provenance indexes; mark missing new-version results as pending validation instead of substituting old data.
- Record `final_selection` (version, confirmation source, scope) and `release_confirmation` (whether the complete parameter snapshot is confirmed) in `design.json`. Selecting a candidate does not confirm every parameter, result, or fabrication revision.

## 7. What Git (GitHub) tracks

The GitHub copy holds what is needed to rebuild and re-simulate the current version; everything else stays in the local working tree.

| Tracked | Local only (`.gitignore`) |
|---|---|
| `pattern design/`, `code/`, `skill/`, `note/`, `templates/`, READMEs | generated GDS: `Design*/gds/**/*.gds` (the layout code regenerates them) |
| per-Design `layout/` (GDS-generating code), `simulation/` (scripts that rerun the simulations) | run folders `Design*/results/*/*/`, except the records listed on the left |
| `design.json` (adopted / candidate results, `fab_junction_targets` EJ / EL), `parameters/*.json` | in kept runs: `code/`, `inputs/`, `raw/`, `logs/` snapshots |
| records of runs referenced by `design.json`: `run.json`, `summaries/`, `figures/`, top-level `*.json` / `*.md` / `*.txt` / `*.png`, and any file `design.json` points to | runs that `design.json` does not reference, or lists under historical / failed / removed keys |
| current `fig/` standard set, `notes/`, `Fabrication output/Template/` marker references | dated fabrication deliveries, AEDT projects, CSV data |

- After changing `adopted_results` or a candidate in `design.json`, run `python -I -X utf8 -B code/general/sync_tracked_results.py` (dry run), then `--apply`: it force-tracks the new records and untracks ignored files. It changes only the Git index; local files are never deleted.
- Links in READMEs and notes to GDS files, raw data or unreferenced runs resolve only in a local checkout.
- Untracking does not shrink Git history; older commits still contain the files removed from tracking.
