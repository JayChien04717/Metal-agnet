# pattern design (shared QComponents)

This directory contains components (QComponent) and port positions shared by Design1–6. Chip layouts remain in each Design's own
`layout/`（`ftf_layout.py`、`purcell_layout.py`、`design5_layout.py`、`chip6_layout.py`…），
Those layouts add this directory to sys.path.

| Category | File | Class |
|---|---|---|
| Port | `pcbport.py` | `pcb_ports` (PCB port positions) |
| Launch pad | `bondpad.py` | `BondPad` / `add_bond_pad` ([Design6 shared parameters](launch_pad.md)) |
| Fluxonium | `my_fluxonium.py` | `MyFluxonium` (default flux line 80 / 11 µm, Design3; Design1 sets 50 / 15 µm in `FX_OPTIONS`) |
| | `fluxonium.py` | `FluxoniumArrayJJ` (all fluxonium simulations; same loop as `MyFluxonium._loop`: JJ is a short segment in the center lead at the qubit center, each U side is a 100 nH lumped L, and the U ends are metal; Q3D uses `flux_as_ground='True'`) |
| | `fluxonium_withJJ.py` | `FluxoniumPocket` (Design2; simulation uses `u_array='True'` for the same loop) |
| transmon / coupler | `my_tmon.py` | `MyTmonCoupler` |
| | `tmon_c.py` | Legacy `Tmon_c` and optional `Tmon_c_B2` using the Design6 short L line |
| | `fluxline.py` | Shared flux-bias-line geometry and Design6 B2 parameters: L bar 45 µm, CPW 10/6 µm, separation 9 µm, x offset +12 µm |
| Capacitor | `hammer_cap.py` | `HammerCap` (C_j) |
| | `enclosing_cap.py` | `EnclosingCap`（hugger，C_in / C_p） |
| | `interdigital_cap.py` | `InterdigitalCap` |
| Other | `tap_junction.py`, `res_end.py`, `test_junction.py` | T junction, open/short termination, and test junction |
| GDS tools (not QComponents) | `code/general/gds_tools.py` | `export_gds` (junction-free), `plot_gds`, `junction_window` (shared by Design1–3 with consistent figures; see `skill/gds_export.md`, `skill/layout_figures.md`) |
| Markers | `xmark.py` | `XMark` (400 µm etched square `window_size`, central 250 × 10 µm metal cross; etched `label` to the right of the window and below the cross, `label_gap` 30 µm, compact character spacing); `add_corner_marks()` (one per corner, centers 0.35 mm from chip edges, die ID at upper left) |
