# Q3D check: readout arc ↔ coupler C-arc parasitic coupling (r6 circular design language)

Date: 2026-11-27. Requested by the user as a small Stage 2 check of the `r6` circular layout.

## Model and settings

- Script: [`simulation/qubit/q3d_ro_coupler_crosstalk.py`](../simulation/qubit/q3d_ro_coupler_crosstalk.py).
- Geometry (sim): D5-like qubit `Q0` at the origin (readout arc NE, JJ SW, XY not drawn), four C-arc couplers `TC_E/N/W/S` around it, four neighbour qubits (no readout) loading the far arcs. Flux lines drawn as ground; JJ landing metal omitted.
- The readout pin is opened by a rotated endcap cut into the `Q0` pocket. qiskit-metal's Q3D `open_terminations` draws axis-aligned endcaps only, which leave a `45°` stem shorted to ground. That is why run `…184042…` failed and was recorded as failed.
- Q3D: `max_passes = 8`, `min_passes = 2`, `percent_error = 0.5`, `Medium` order, `5 GHz`. Each run took about 2 min. Pass count is an execution setting; convergence was not checked separately (status: single run).
- Readout node = arc pad only, so `g_qr` / `g_cr` are pad-level estimates (λ/4 `C_r = π/(4 ω_r Z0) = 349 fF`, `f_r = 7.17 GHz`, `f_q = 4.65 GHz`, `f_c = 6.2 GHz`). The ratio `g_cr / g_qr` is the robust number.
- Sanity checks (`skill/cap_matrix_analysis.md` §7) passed in every completed run: symmetric matrix, positive diagonal, non-positive off-diagonal, positive row sums.

## Results

| Run (`results/qubit/`) | C-arc span | RO arc | Spokes | C(RO,Q0) fF | C(RO, adjacent C-arc) fF | g_qr MHz | g_cr (adjacent) MHz | g_cr / g_qr | g_qc MHz | Ec Q0 / TC MHz |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|
| `20261127_184405_…_r6_0666cb1a` (r6 as drawn) | 36° | 26° | no | 3.44 | 0.76 | 56.6 | 7.87 | **13.9 %** | 109–111 | 219 / 262 |
| `20261127_184703_…_c28_ro26_bfe0410f` | 28° | 26° | no | 3.51 | 0.49 | 57.8 | 5.10 | **8.8 %** | 94–95 | 219 / 265 |
| `20261127_184950_…_c28_ro20_b3364fb3` | 28° | 20° | no | 2.95 | 0.36 | 48.5 | 3.75 | **7.7 %** | 94–95 | 219 / 265 |
| `20261127_185909_…_c28_ro26_spokes_340095fd` | 28° | 26° | yes (2 × 10 µm) | 3.12 | 0.36 | 51.0 | 3.77 | **7.4 %** | 90–94 | 217 / 265 |

- Only the two couplers whose arcs flank the readout arc couple noticeably (`C ≈ 0.4–0.8 fF`). The opposite couplers are at `0.02 fF` (`g_cr/g_qr ≈ 0.4 %`), and the far-side arcs (`pad_r`) at `< 0.001 fF`.
- The coupling runs across the shared etched pocket (line of sight between arc ends). Narrowing the arcs is the main lever. Grounded shield spokes add only a little (`8.8 % → 7.4 %`, partly by also lowering `g_qr`), because the field also links through the substrate and vacuum above the plane.

## Impact estimate (r6 as drawn, `g_cr = 7.9 MHz`, `Δ = f_c − f_r ≈ −0.97 GHz`, `α_c ≈ −0.26 GHz`)

- Coupler-state-dependent resonator shift `χ_rc ≈ g² α/(Δ(Δ+α)) ≈ −14 kHz` (≈ 2 % of the readout `χ ≈ 0.75 MHz`). With `C-arc 28°`: ≈ −3 kHz.
- Coupler Purcell decay through the (unfiltered, `κ_r = 2 MHz`) resonator: `(g/Δ)² κ_r ≈ 133 Hz` → `T1,c ≈ 1.2 ms`. With `C-arc 28°`: ≈ 5 ms. The Purcell filter improves both.
- Readout photons Stark-shift the coupler by `2 χ_rc n̄`. At `n̄ ≈ 20` that is ≈ 0.5 MHz for `r6` and ≈ 0.1 MHz for `C-arc 28°`. This slightly moves the ZZ-zero idle point during readout.
- Conclusion: not catastrophic, but `r6` (`13.9 %`) leaves little margin. `C-arc 28°` halves the parasitic and brings `g_qc` closer to its `70–85 MHz` target.

## Other findings from the same matrices (targets for Stage 2 tuning)

- `g_qc` (qubit–coupler) is larger than the target (`110 MHz` at 36°, `94 MHz` at 28°, target `70–85 MHz`).
- `Ec(Q0) ≈ 219 MHz` (target `190–200 MHz`). `Ec(TC) ≈ 262–265 MHz` (target `180–200 MHz`). Both islands need more total capacitance.
- The qubit mode is not distorted asymmetrically: the four `g_qc` agree within 2 %, and the readout-side couplers differ only by the small RO loading.

**Adopted (user decision, 2026-11-27):** `C-arc 28°`, readout arc `26°`, no shield spokes. Applied to layout `r6` (`COUPLER_OPT.c_span_deg = 28`). The matching evidence is run `20261127_184703_…_c28_ro26_bfe0410f`. The `ro_shield_spokes` qubit option remains available, default off.
