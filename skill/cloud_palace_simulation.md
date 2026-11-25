# Skill: Cloud EM Simulation with AWS Palace (`Qiskit Metal` → `SQDMetal` → `Gmsh` → `Palace`)

Design-agnostic operating rules for running electrostatic (capacitance), eigenmode, and frequency-driven (S-parameter) simulations in Google Colab Pro (`CPU + High-RAM`) with `SQDMetal` + `Gmsh` + `AWS Palace`.

---

## 1. Cloud Environment & Imports

* **Runtime:** Google Colab Pro **CPU + High-RAM** (8 vCPU / ~51 GB RAM). Multiple notebooks can run in parallel on independent VMs; write outputs to timestamped subfolders under `runs/<RUN_NAME>_<timestamp>/` so concurrent runs never collide.
* **Bootstrap:** Mount Drive once, then run `_env/setup_colab_env.sh` (restores precompiled `/opt/palace/bin/palace` and `/content/sqd_env` Python 3.11 environment in ~1–2 min).
* **No AEDT / No markers in cloud:**
  - Guard `code/general/project_paths.py` `_install_aedt_session()` so missing `code/aedt/aedt_session.py` returns silently instead of raising `FileNotFoundError`.
  - Stub GDS marker modules before importing chip layouts that reference them (cloud EM runs never need wafer alignment marks or `gdspy`):
    ```python
    import sys, types
    sys.modules.setdefault('xmark', types.SimpleNamespace(add_corner_marks=lambda *a, **kw: None))
    ```
  - Import shared `QComponent`s directly from `pattern design/` and parameters from JSON rather than legacy wrapper scripts.
* **Error visibility in notebooks:** Invoke Python scripts with `subprocess.run(..., capture_output=True, text=True)` and write stdout/stderr to a `.log` file on Google Drive before checking `returncode != 0`, so failures are always preserved even if Colab output frames disconnect.

---

## 2. Geometry Preparation & SQDMetal / Gmsh Pitfalls

### 2.1 Gmsh OCC `"Could not create line"` on curved Shapely polygons
* **Symptom:** `gmsh.model.occ.add_line(points[-1], points[0])` raises `Exception: Could not create line` inside `GMSH_Geometry_Builder._draw_polygon_in_GMSH_from_coords`.
* **Cause:** Shapely `buffer()` / `unary_union()` rounded corners leave near-coincident consecutive or closing vertices (`< 1e-6 mm`), which Gmsh OCC rejects as zero-length edges.
* **Fix (mandatory for layouts with rounded brackets/fillets):**
  1. Simplify complex polygons with `geom.simplify(5e-4, preserve_topology=True)` (0.5 µm tolerance).
  2. Monkey-patch `GMSH_Geometry_Builder._draw_polygon_in_GMSH_from_coords` to deduplicate consecutive and wrap-around points closer than `1e-6 mm` before calling the original method.

### 2.2 Open terminations in capacitance runs (`add_ground_plane()` shorting trap)
* **Symptom:** A 2-electrode component (e.g., an enclosing or interdigitated capacitor coupon) collapses to a single `1×1` capacitance entry in `terminal-C.csv`.
* **Cause:** Unlike Q3D's `open_terminations=[(comp, pin)]`, SQDMetal's `add_ground_plane()` merges any metallic trace touching the ground plane into the ground net.
* **Fix:** Explicitly place `OpenToGround` components on every open pin (or shorten leads so they terminate strictly inside the etch pocket with a nonzero gap to ground). Always verify that `terminal-C.csv` has `N_conductors + 1` rows/columns (ground + all intended signal islands) before trusting the matrix.

### 2.3 Substrate permittivity ($\epsilon_r$) & thickness matching
* **Critical:** SQDMetal's built-in `'silicon'` defaults to **$\epsilon_r = 11.45$**, whereas standard Ansys HFSS/Q3D silicon uses **$\epsilon_r = 11.9$** (or `11.5` in specific presets).
* For CPW structures ($\epsilon_{\text{eff}} \approx (\epsilon_r + 1)/2$), $\epsilon_r = 11.45$ vs $11.9$ shifts:
  - Capacitances by **$-3.5\%$** ($E_C$ by **$+3.6\%$**),
  - Resonator frequencies by **$+1.8\%$**.
* Always record `Permittivity` from Palace's generated `<run>.json` and match substrate thickness (`size_z`) and $\epsilon_r$ explicitly when comparing against local AEDT runs.

---

## 3. Meshing, Solver Order, and Speed Guidelines

| Simulation type | Recommended domain | `min_size` rule | `solver_order` | Runtime guidance (Colab 8 vCPU / 51 GB) |
|---|---|---|---|---|
| **Capacitance (Electrostatic)** | Crop chip to active core + ≥ 0.4 mm ground buffer | **≤ 0.5× smallest critical gap** (e.g. `2–4 µm` for 4 µm gaps; `10–12 µm` for ≥ 10 µm gaps) | `2` (production), `1` (draft) | 20 s – 3 min |
| **Eigenmode** | Single resonator + coupler + qubit (crop to active region + ≥ 0.5 mm buffer) | `10–15 µm` on CPW/electrodes, `max_size = 120–160 µm` | `2` (or `1` for quick scan) | 5 – 15 min |
| **Frequency-Driven S21** | Single channel + `LaunchpadWirebond` CPW ports (`create_port_CPW_on_Launcher`) | `10–15 µm` on CPW/resonator | **Never** run 30+ uniform `Order: 2` points (~20 min/pt!) | Use **Adaptive Frequency Sweep (`AdaptiveTol: 1e-3`)** or `Order: 1` coarse scan first |

### Driven S21 optimization (avoiding OOM, disk exhaustion, and 10+ hour sweeps)
* **51 GB RAM ceiling (`Order = 2` OOM trap):** Combining `solver_order = 2` with a fine mesh (`min_size ≤ 10–14 µm`) and `4 MPI ranks` (`SuperLU`) exceeds Colab Pro's 51 GB RAM limit and aborts with `mpirun ... exited on signal 9 (Killed)` before writing `port-S.csv`.
* **CRITICAL `"SaveStep": 1` disk-exhaustion trap:** SQDMetal translates `solns_to_save: 1` into `"SaveStep": 1` in `Solver.Driven`, which tells Palace to write a full 3D ParaView `.vtu` dump **at every frequency point** (~200 MB/pt × 1,001 pts = **~200 GB**, crashing Colab disk and bloating Drive sync!). Always patch `j['Solver']['Driven']['SaveStep'] = 0` (or `500` to keep only 2–3 snapshots) before calling `s21_sim.run()`, and clear old `/content/palace_runs/s21_driven/` folders between runs.
* **Recommended Colab S21 recipe:**
  1. Use **`solver_order = 1`** with fine mesh (`8–10 µm` on tight capacitor/CPW gaps) and **`num_cpus = 2, num_threads = 4`** (halves MPI-duplicated `SuperLU` memory).
  2. Always enable Palace's **PROM Adaptive Frequency Sweep** (`"AdaptiveTol": 1e-3`, `"AdaptiveMaxSamples": 20–50`, **`"SaveStep": 0` or `500`**) in `Solver.Driven` so Palace solves only ~10–25 full-order snapshots and interpolates 0.5–2 MHz fine curves in **15–35 minutes** instead of 10+ hours.
  3. For multi-channel chips (e.g. 6 channels / 12 poles), set `AdaptiveMaxSamples ≥ 8 × N_channels` or split into ~0.5–0.7 GHz sub-bands so narrow-linewidth resonances are not skipped.
  4. Remember that planar CPW layouts without ground wire bonds / air bridges can excite parasitic slotline modes on large chips.

---

## 4. Plotting & Post-Processing (`common/plot_utils.py`)

* **Unit check:** Palace ParaView `.pvd`/`.vtu` coordinates use the model length scale `L0 = 0.001` (meters → **mm** in ParaView output, or meters depending on SQDMetal version). Verify coordinate range before multiplying by `1e3` so axis labels (`mm` vs `µm`) match the plotted numbers.
* **Capacitance post-processing:** Palace outputs a full `(N+1) × (N+1)` Maxwell matrix in `outputFiles/terminal-C.csv` (including the ground plane as conductor 0, in Farads). Drop the ground row/column and apply the standard inverse-matrix procedure (`skill/cap_matrix_analysis.md`) to extract $E_C$, $J$, $g$, and $\beta$.
* **Eigenmode verification:** Always inspect the `log10(|E|)` field plot (`time_point = m - 1`) to confirm which eigenmode is the physical CPW resonator mode vs a box/substrate mode before quoting $f_r$.
