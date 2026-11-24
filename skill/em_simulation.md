# Skill: EM simulation environment, solver defaults, and wire bonds

Shared execution rules for Q3D and HFSS runs across all Designs. Design-specific solver overrides (for example, fewer passes on a large multi-channel chip) belong in that Design's `README.md` / `notes/` and must be recorded in `run.json`.

- Capacitance matrix extraction: [cap_matrix_analysis.md](cap_matrix_analysis.md)
- Mutual inductance extraction: [mutual_inductance_analysis.md](mutual_inductance_analysis.md)
- Purcell filter + readout channel validation: [purcell_resonator_design.md](purcell_resonator_design.md)
- Run recording and provenance: [parameter_provenance.md](parameter_provenance.md)

---

## 1. Python and AEDT environment (do NOT search for Python every time)

* **Layout, GDS, and Q3D/HFSS drivers:** conda env `qmetal` (`qiskit_metal`, `gdspy`, `geopandas`, `shapely`, `pyEPR`).
  - Lab PC: `C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe`
  - Laptop (user `zx123`): `C:\Users\zx123\miniforge3\envs\qmetal\python.exe`
  - AEDT 2021 R2 on both. Do not use `.venv` or system Python (they lack `geopandas` / `gdspy`).
* **Quantum analysis scripts:** use the environment documented by the target Design (e.g. `scqenv`).
* Run scripts from the repository root or the Design's `layout/` / `simulation/` folder per the script's entry point:
  `& "C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe" -I <script>.py`.
* `qiskit_metal` prints a harmless `Fake_Junctions.GDS` warning on stderr (PowerShell may show a `NativeCommandError` header); ignore it and read the output below it.

---

## 2. Default pass counts and initial mesh

| Solver | Max passes | Initial mesh (length based) |
|---|---|---|
| **Q3D** (capacitance matrix) | **8** | Q3D auto |
| **HFSS eigenmode** | **6** | **10 µm** on every CPW line and capacitor electrode |
| **HFSS driven S21** | **4** | **Explicit length-based mesh required** (default **10 µm** on every CPW line and capacitor electrode) |

* **Why:** seeding an explicit length-based mesh on the metal from pass 1 avoids slow pass-by-pass refinement and keeps memory bounded.
* **Explicit mesh length required for HFSS driven:** always specify the initial mesh length in the script and record the mesh-length value and target faces/objects in `run.json`.
* **Deviations:** deviate only deliberately (e.g., low-frequency mutual-inductance Z-parameter runs, high-precision κ checks, or large multi-channel full chips that exceed RAM). Document Design-specific overrides in that Design's `README.md` / `notes/` and in `run.json`.
* **Full 10 × 10 mm multi-channel chips:** a 10 µm seed mesh across many channels does not fit in 24–32 GB RAM. Tune per channel or channel pair; run the full chip only with an explicitly documented coarser override.

### HFSS driven acceptance

* On broadband or well-damped driven runs, a final adaptive maximum magnitude ΔS ≤ **0.05–0.10** is acceptable.
* On narrow lossless resonances (e.g. Purcell channels), adaptive ΔS typically does not settle pass-to-pass (it jumps 0.08–1.8 on a −60 dB null). Keep the fixed pass budget (4 passes with explicit length mesh) and verify the extracted physical quantities separately.
* Pass counts are execution settings, not convergence evidence (`parameter_provenance.md`). Always record the final ΔS / ΔF, the configured solver target, and whether the run met that target.

---

## 3. Ground wire bonds in every HFSS / Q3D run

* Enable `hfss_wire_bonds = True` on every etched path in the path QGeometry table and set
  `sim.renderer._options.update(wb_threshold='350um', wb_offset='0um', wb_size=2)` before `run`.
  (The default `400um` skips most 400 µm serpentine turns; `350um` places a bond in each serpentine bottom and every ~300–350 µm on long runs.)
* Short etched pieces below the threshold (coupler stubs, tee arms) need extra PEC straps drawn after `render_design`.
* Use the same bonds in Q3D, eigenmode, and driven runs. Q3D treats bond nets as ground; if a bond net fails to touch ground on both sides, Q3D fails (run Q3D without that bond—it barely affects series C).
* Plot bond positions before a long run.

---

## 4. HFSS practical lessons

* **Mesh Link from eigenmode:** where possible, build an eigenmode design with the same geometry and 50 Ω RLC port sheets, then link a driven design in the same project (`Mesh Link -> <eig design> : Setup : LastAdaptive`, 1 pass). AEDT can refuse a Mesh Link silently (`ImportMesh` stays `False`)—always verify it.
* **Eigenmode basis order:** use first-order basis; second order exhausts 32 GB RAM at pass 1.
* **Locate resonances first:** an Interpolating sweep can miss a narrow dip. Run a wide sweep (e.g. 5.8–7.6 GHz, 1 MHz step) first, then a narrow sweep (≤ 0.25 MHz) around each dip. Discrete 5 MHz steps miss filter dips.
* **Far-detuned resonators (kHz linewidth):** invisible in a driven sweep; locate them with eigenmode.
* **Scratch disk space:** keep ≥ 80 GB free on the HFSS temp drive (`%TEMP%\hfss_*.pjt` reaches 20–45 GB per session). A full disk, a disconnected net, or OOM all raise COM error `-2147352567` at `Analyze`—read the AEDT message window before guessing.
* **AEDT session cleanup:** an AEDT that a run launches is closed when the run ends; an AEDT that was already open is left alone (`code/aedt/aedt_session.py`, installed automatically when `project_paths` is imported, so scripts need no extra code). At exit the run saves and closes every open project and quits AEDT; a detached watchdog force-kills that AEDT and its solver children if it is still running 90 s after Python exits, which also covers a run that was hard-killed or timed out. For interactive work where AEDT should stay open between runs, open AEDT yourself first or set `AEDT_KEEP_OPEN=1`. `python code/aedt/aedt_session.py --list` lists running AEDT processes; the watchdog log is `%TEMP%edt_session_watchdog.log`. A force-killed AEDT can leave `<project>.aedt.lock`; delete it only when no AEDT is running. Agents must not leave idle AEDT sessions behind.
* **Shared project & side files:** all runs of a Design share one AEDT project (`save_aedt.use_project`). Never leave solver side files (such as `hfss_eig_f_convergence.csv`) in the repository root; `project_paths.configure()` and `run_record.create_run()` redirect them into the active run's `raw/` folder.
