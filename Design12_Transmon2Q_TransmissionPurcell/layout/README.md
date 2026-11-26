# Design12 layout

- `build.py` is the single entry point: `build_design(variant='sim', overrides=None, channels=True, flux_as_ground=False)`.
  - Running it as a script prints the readout length report and writes `gds/Design12_Transmon2Q_TransmissionPurcell_2026-11-26_r1_sim.gds`.
  - Override keys: `pf`, `idc`, `rr`, `xmon`, `xmon_Q1` / `xmon_Q2`, `ro_claw`, `ro_claw_Q1` / `ro_claw_Q2`, `coupler`, `qx_mm`.
- **Layout constants:** as in Design11, the knobs are module-level dicts in `build.py`, not `placement.json` / `routes.json`. `parameters/geometry.json` mirrors the `r1` values.
- **Lengths:** filter and resonator lengths are solved in closed form in `build.py` (fixed bar count; meander width or bus position adjusted). The drawn filleted lengths are checked against the targets.
- **Components:** `components/` holds Design-specific QComponents (see its README).
