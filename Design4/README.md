# Design4 — Purcell filter + readout resonator


> **Directory reorganization (2026-11-15)**: [Layout code](layout/) · [Simulation entry points](simulation/README.md) · [Results index](results/README.md) · [Adopted version](design.json). Solvers were not rerun. Moving historical values does not make them final confirmed values. Create a new run for new executions and confirm final parameters before commit/push.


The current layout uses the [Design6 launch pad standard](../skill/launch_pad.md). Saved solver results correspond to the launch geometry used in their original runs.

## Importable PurcellResonator

`PurcellResonator` packages the Purcell filter, readout resonator, and feedline using 10/6 µm CPW.
The hammer electrodes and leads on both sides have been removed, retaining the open stubs and existing HFSS ground-strap configuration.
Connect `feed_in` and `feed_out` to routes or ports in the host layout.

```python
import sys
from pathlib import Path

from qiskit_metal import designs

# Run from the SQL HFSS directory, or use the absolute path to Design4/layout.
sys.path.insert(0, str(Path("Design4/layout").resolve()))
from purcell_resonator import PurcellResonator

design = designs.DesignPlanar()
design.chips.main.size.size_z = "-500um"
pair = PurcellResonator(design, "readout_channel", options={
    "cpw_spacing": "400um",
    "resonator_trim": "0um",
})
```

- `cpw_spacing` is the center-to-center spacing of the two parallel CPWs. The default is 400 µm; values below 250 µm produce a warning. In the no-hammer HFSS eigenmode sweep, the minimum sampled mode splitting was 2.56 MHz at 250 µm, 1.34 MHz at 300 µm, and 4.10 MHz at 350 µm. The 400 µm reference gave 5.06 MHz. The warning threshold uses the shortest sampled spacing meeting the approximately 5 MHz target. Sweep data: `results/resonator/legacy/2026-11-11/no_hammer_spacing/`.
- `filter_length` and `resonator_length` set the nominal route lengths; `resonator_trim` moves only the resonator end. Default lengths come from the spacing-comparison model. Validate frequencies for the actual host chip and boundary conditions.
- `ground_wire_bonds` defaults to `True`, marking paths for automatic HFSS ground straps. The host design sets chip dimensions, materials, and physical fabrication bridges.
- HFSS spacing-comparison results are in `no_hammer_eigen` and `no_hammer_400_eigen` under `results/resonator/legacy/2026-11-10/purcell_6500_paper_hammer/`.

## Connect the ground around the hammer (2026-11-11)

The etch of `stub_f` + hammer pocket + `stub_r` forms a continuous cut, dividing the ground between the filter and resonator into upper and lower regions.
HFSS consequently shows an additional ground-current coupling path: at 190 µm, the splitting remains 122 MHz without hammer electrodes if the pocket remains continuous.
This path is not capacitive and is absent from Q3D; fitting assigns its contribution to C_j (10–11 fF, versus 2.4–2.6 fF from Q3D).
This is the purpose of the air bridges at the T-junctions on both sides of Wallraff's hammer (C_AB in Fig. 2.29).

- `purcell_6500.hammer_ground_straps(design)`: one ground strap at the center of each stub (60 µm span, 15 µm wire diameter); landing points are checked automatically to avoid etched regions.
- `purcell_6500.draw_hammer_straps(renderer, design)`: call after `renderer.render_design(...)`; draws PEC bond wires named `g_wb`, as in qiskit-metal auto_wirebonds. Q3D still treats them as ground.
  qiskit-metal auto_wirebonds places straps only on straight segments longer than 400 µm, so these two stubs previously had none.
- Use `build(..., spacing_mm=0.25)` to change the spacing directly. Illustration: `fig/purcell_6500_hammer_straps.png`.

**Ground straps are required on every CPW segment**, including hammer stubs, T-junctions, hugger stubs, and short meander segments.
HFSS eigenmode at 250 µm spacing and 9 µm gap: J = 68.2 MHz with qiskit auto bonds alone; 22.8 MHz after adding 2 hammer straps;
20.4 MHz with straps on every CPW segment (62 straps, `ground_straps(design, 'cpw')`), corresponding to C_j = 3.04 fF versus 2.40 fF from Q3D.
Script: `simulation/resonator/hfss_hammer_straps.py --straps cpw`. Details and data: `notes/purcell_design4.md`.

## Verification records for the former Design6


(Starting 2026-11-09. The former Design6, with Cai 2024 channels and c6_* runs, was deleted; recover it from Git history if needed.)

## Questions to answer

1. **Purcell ↔ feedline**: do ω_p and κ_p predicted by the circuit model using Q3D C_p agree with HFSS S21?
2. **Resonator ↔ Purcell**: does J predicted using Q3D C_j agree with HFSS S21?

Conversely, how far do C_p and C_j fitted to HFSS S21 differ from their Q3D values?

## Design rationale

- **No C_in**: the feedline connects the left and right pads directly and is matched at both ends. S21 = 2 / (2 + Z0·Y_branch), apart from a phase factor.
  Standing waves do not modify κ_p (with C_in, κ_p varies from 0.9–93 MHz with position; see Design5 README 9b).
- **Verify one coupling at a time**:
  - `filter`: C_p + Purcell filter only. Remove the resonator route but retain the C_j hammer with its opposite electrode floating, preserving the full-model filter loading.
    S21 has **one** notch, giving ω_p and κ_p directly; these depend only on C_p and the filter length.
  - `full`: complete channel. C_p was verified in the preceding step; the remaining unknown is C_j → J.
- One channel (6.7 GHz, using Design5 chip6 ch3 lengths), with hairpin routing.

## Files

```
Design6/
├─ README.md
├─ layout/verify_layout.py   layout: VARIANT filter / full (sweep capacitance using CP_PLATE and CJ_GAP environment variables)
├─ simulation/resonator/verify_model.py        circuit model without C_in and predictions directly from capacitances
├─ gds/  verify_filter.gds, verify_full.gds
└─ fig/  verify_filter(_zoom).png, verify_full(_zoom).png
```

## Procedure

| Step | Contents | Comparison |
|---|---|---|
| 1. Q3D | Solve the entire `verify_full` channel | C_p, C_j, and ground capacitances for this geometry without C_in |
| 2. HFSS `filter` | S21/S11 | Q3D-capacitance model versus HFSS: ω_p, κ_p; fit ΔL_filter and the C_p scale → fitted C_p versus Q3D |
| 3. HFSS `full` | S21/S11 | Fix C_p to step 2; fit ΔL_res and the C_j scale → fitted C_j versus Q3D, and Q3D J versus fitted J |
| 4. Optional sweep | 3 C_j gaps and 3 C_p plate sizes | Check J ∝ C_j and κ_p ∝ C_p² |

**HFSS workflow for the Purcell-filter project**: mesh with eigenmode first (first order, initial mesh of 10 µm on lines / 5 µm on capacitors,
Delta F 0.1 %). Use Mesh Link in the driven solve to reuse that mesh, solve only 1 pass without refinement, then sweep
(`simulation/resonator/hfss6_meshlink.py`). A filter-only HFSS solve is unnecessary: without C_in, the full-channel S21 can simultaneously fit
ΔL_filter、ΔL_res、C_p、C_j。

Each run outputs only capacitances (Q3D and fitted), HFSS S21/S11, and one comparison figure (raw HFSS / Q3D-capacitance model / fitted model).
Store in `results/resonator/legacy/<YYYY-MM-DD>/<HHMM>_<tag>/` without overwriting.

## Current predictions (Design5 full-channel Q3D: C_p 45.91 fF, C_j 3.43 fF; step 1 replaces these with values for this geometry)

| | ω_p (GHz) | ω_r (GHz) | J (MHz) | κ_p (MHz) |
|---|---|---|---|---|
| filter only | 6.7202 | — | — | 37.4 |
| full | 6.7190 | 6.7000 | 19.3 | 37.3 |


## Design notes

Design records, including tuning history, Q3D/HFSS data, and implementation: [`notes/purcell_design4.md`](notes/purcell_design4.md).
