# Flux line → fluxonium loop mutual inductance (HFSS driven modal)

Status: simulated for the current FX1 geometry (see *Results*): M = 0.44 pH.
Geometry `pattern design/fluxonium.py` (`FluxoniumArrayJJ`), script
`simulation/hfss_mutual_fx.py`.

## Geometry

`FluxoniumArrayJJ` = `MyFluxonium` (same pads, pocket, U-loop flux line, readout and
charge line) with

- the hairpin array at its real size: outer 55 × 10 µm, 2 µm wires. As in
  `FluxoniumPocket` (`qiskit_metal/qlibrary/qubits/fluxonium_v2.py`), the inductive
  parts are junction-table lines with an inductance. The two sides of the U,
  `array_top` / `array_bot`, are lumped L (2 µm wide, `l_inductance / 2` = 100 nH each), joined by
  a metal closing end `array_end`. Port-only JJ: the L values do not change M.
  Before: two lumped arms of `l_inductance / 2` joined by a metal end;
- one 2 µm lead line from the pads (stems);
- the main JJ: the 2 µm gap in the middle of the lead line, bridged by `rect_jj` (`L_j`, or
  Port 2), at the pocket centre. Before 2026-11-12 the JJ sat on a separate branch 6 µm right
  of the leads; the array and flux tip have moved 6 µm with the new loop, so M has to be re-run.

With Port 2 as port only (Z parameters, I2 = 0) no current flows around the loop, so
the lumped L values do not change M; they matter for eigenmode or other driven runs. Film penetration / kinetic inductance of the 2 µm wires is not modelled.

Layout plots: `fig/mutual/fluxonium_arrayjj.png`, `fig/mutual/fluxonium_arrayjj_zoom.png`
(`python fluxonium.py`).

## Ports

- Port 1: the `flux_line` pin (pocket wall, top leg of the U; the bottom leg is shorted
  to ground at the wall). `port_list=[('FX1', 'flux_line', 50)]`; the renderer cuts a
  `cpw_gap` endcap and puts the lumped port across it.
- Port 2: the JJ gap. Use `jj_to_port=[('FX1', 'rect_jj', 50, False)]` (port only):
  the port voltage is then the open-loop EMF of the whole loop. `True` draws L_j in
  parallel with the port; that closes the loop through L_j + L_array and the port sees
  only EMF · L_j / (L_j + L_array) (≈ 0.15 for 34 nH / 200 nH), so Im Z21 / ω would
  underestimate M by that factor.

## Extraction

V2 = jω M I1  →  **M = Im Z21 / (2π f)**, HFSS report `im(Z(2,1))/(2*pi*Freq)`.

Valid in the quasi-static range, far below any structure resonance. Start with
10, 30, 100 MHz and check that M is flat in f; too low a frequency makes the pH-level
signal hard to resolve numerically.

The flux per unit bias current is then Φ/I = M, i.e. I_Φ0 = Φ0 / M.

## Script: `simulation/hfss_mutual_fx.py`

- Renders one `FluxoniumArrayJJ` with the FX1 options (FX_OPTIONS + FX_PER_QUBIT) on its
  own; ground plane = bounding box + 0.2 mm buffer (`box_plus_buffer=True`), 525 µm Si.
  The flux line here returns through the local ground only; the chip routing and
  launch pad are not included.
- Readout and charge-line pins are open-terminated (endcap gap).
- Setup: adaptive at 100 MHz, max ΔS 0.0002, 12–15 passes (≥ 12 forced), 3 converged;
  discrete sweep 10–100 MHz, 10 points (~2.6 min). ΔS is a weak criterion at these
  frequencies: max ΔS 0.002 stopped after 4 passes (7k elements) with M 7 % low.
- Output in `simulation/results/runs/`: `<tag>_Z.csv` (per frequency: M, L_fl = Im Z11/ω,
  L_loop = Im Z22/ω, Re/Im Z21, Z21–Z12 mismatch), `<tag>_summary.json`,
  `<tag>_convergence.csv`.
- Checks: M flat over the sweep (spread printed); Z21 ≈ Z12; L_loop ≈ l_inductance
  (200 nH) + geometric loop inductance.
- The JJ port is the one with the larger Im Zii (~200 nH vs ~80 pH). The HFSS
  `GetExcitations` order does not match the exported Z-matrix order, so it is not used.

## Results (2026-11-03, FX1 options, fl_sep 20 µm, l_inductance 200 nH, L_j 34.38 nH)

| run (`simulation/results/runs/`) | passes / elements | M @ 10 MHz | L_fl | L_loop |
|---|---|---|---|---|
| `mutual_fx` (old setup, ΔS 0.002) | 4 / 7.3k | 0.409 pH | 62.2 pH | 201.1 nH |
| `mutual_fx_conv` | 12 / 63k | 0.4414 pH | 79.5 pH | 201.7 nH |
| `mutual_fx_p15` | 15 / 139k | 0.4415 pH | 79.6 pH | 201.7 nH |

- **M = 0.44 pH → I_Φ0 = Φ0 / M ≈ 4.7 mA** (bias current for one flux quantum).
- Z21 = Z12 to < 1e-7; Re Z21 ≈ 0. L_loop ≈ 200 nH array + 1.7 nH geometric/lead.
- M rises by 1.6 % from 10 to 100 MHz (∝ f², onset of the loop / pad resonance), so
  the 10 MHz value is quoted as the quasi-static M.
- Model limits: PEC sheets (no kinetic inductance / film penetration of the 2 µm
  wires), flux current returns only through the local ground (no chip routing).
- Cross-check: free-space Biot–Savart of thin filaments (U-loop centre lines vs the
  inner hairpin + JJ-branch area, no ground plane) gives 0.44 pH.

## Flux-line tip distance (2026-11-03, current design)

Same structure (U-loop from the side wall, tip facing the array's closed end); leads now
reach the arms' inner edge (full contact with the lumped-L sheets). 12 passes.

| run | tip – array end | M @ 10 MHz | I_Φ0 |
|---|---|---|---|
| `mutual_fx_fl20` | 20 µm | 0.440 pH | 4.7 mA |
| `mutual_fx_fl15` | **15 µm (chip layout)** | **0.575 pH** | **3.6 mA** |

The chip layout uses this geometry (now the `MyFluxonium` defaults): array 55 × 10 µm,
2 µm wires and leads, tip 15 µm. `FluxoniumArrayJJ` takes the same options. Q3D of the chip: `simulation/results/latest/fxarray_fl15`.
A grounded line along the array's long side (~1.2 pH free-space) was drawn but not
adopted (it needs the array outside the pads) and has been removed from `fluxonium.py`.

## Reading the impedance response (`fig/mutual/mutual_fx_impedance.png`)

`python simulation/plot_mutual_z.py` (broadband run `mutual_fx_broadband`, 0.01–6 GHz
interpolating, + the converged 10–100 MHz run `mutual_fx_p15`).

- Im Z21 rises linearly with f (slope 1 on log-log) at low f: Z21 = jωM. M is read where
  Im Z21/ω is flat (plateau, 10–100 MHz).
- Pole at **f_r ≈ 2.65 GHz**: L_loop (≈ 200 nH, port 2 open) resonating with the pad
  capacitance; Im Z22 changes sign there. Above it Im Z21 is negative.
- Between ~0.2 GHz and f_r, Im Z21/ω rises faster than the pure inductive
  M / (1 − f²/f_r²): the extra part is capacitive coupling flux line → pads, which grows
  ∝ f² because at low f the 200 nH array shorts the pads. That is also the +1.6 % from
  10 to 100 MHz; at 10 MHz it is negligible, so M(10 MHz) is the inductive M.
- Re Z21 ≈ 0 (lossless), Z21 = Z12 (reciprocal) — sanity checks.
- The kink near 1 GHz in the broadband curve is an interpolating-sweep artefact
  (the discrete 10–100 MHz run is the reference).

## Comparison: reference FluxoniumPocket (`qiskit_metal/qlibrary/qubits/fluxonium_v2.py`)

Script `simulation/reference/mutual_fluxonium_v2.py`, same `run_design` and setup
(layout plot `fig/mutual/fluxonium_v2_mutual.png`);
default options + `make_fbl=True` (pocket 900 × 550 µm, pads 15 × 100 µm + r 60 µm
circles, 30 µm JJ gap, flux U-loop tip at x = 75–85 µm, legs 50 µm apart).

| | M @ 10 MHz | I_Φ0 | L_fl | L_loop | passes / elements | Biot–Savart (free space) |
|---|---|---|---|---|---|---|
| FluxoniumPocket (`mutual_fluxonium_v2`) | **1.19 pH** | 1.74 mA | 518 pH | 200.2 nH | 12 / 80k | 0.84 pH |
| FluxoniumArrayJJ (`mutual_fx_conv`) | 0.44 pH | 4.7 mA | 79.5 pH | 201.7 nH | 12 / 63k | 0.44 pH |

- The reference loop is ~4× larger (inductor 130 µm long, 10 µm from the pad edge:
  ~1500 µm² vs ~360 µm²), which outweighs its flux line being farther away (57 µm
  vs 20 µm). HFSS is 40 % above the free-space estimate there; likely flux focusing
  by the wide PEC pads next to the 10 µm slot (not checked further).
- Z21 = Z12 to < 1e-9; M rises 2 % from 10 to 100 MHz, as for FluxoniumArrayJJ.
- `fluxonium_v2.py` as it stands: `make_pocket` reads `p.nanowire` (not an option;
  the option is `nanowire_inductor`), so the JJ-array branch (130 µm at x = 18 µm) is
  always drawn — that is the geometry simulated. Its `flux_bias_line` pin is 60 µm
  wide (points from y = 0 to −60 µm), so the script re-adds it 10 µm wide with
  gap = cpw_gap (11.233 µm) at the end of the bottom leg (`FluxoniumPocketSim`).

```
python hfss_mutual_fx.py mytag
python -c "from hfss_mutual_fx import run; run('fl10', fx_options=dict(flux_line_options=dict(fl_sep='10um')))"
```
