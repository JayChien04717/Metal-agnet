# r9: Manhattan-compliant JJs and the rounded-rectangle rail-SQUID coupler

Date: 2026-11-28. Stage 2 small-model checks requested by the user. Scripts and run layout as in
[coupler_ring_hub_flux_2026-11-27.md](coupler_ring_hub_flux_2026-11-27.md).

## 1. Fab constraint: Manhattan JJ orientation

- **The evaporation:** Manhattan JJs come from a two-angle shadow evaporation, first along y and then along x. The overlap of the y-finger and the x-finger is the junction, so the fingers must lie strictly along x and y (`pattern design/ebeam_junction/transmon_junction_manhattan.py`).
- **r6–r8 were not compliant:** they rotated every JJ by −45°, which makes the fingers diagonal.
- **What is allowed (user, 2026-11-28):** 0 / 90 / 180 / 270° on one chip, with no extra per-orientation test coupons. A 90° JJ only swaps which island gets the bottom (oxidised) electrode. Never 45°.
- **r9 orientations:**
  - qubits: vertical JJs at a free diagonal corner; north corners use the 180° pattern;
  - couplers: pattern orientation = coupler orientation (0 / 90 / 180 / 270°).
- **HFSS rendering:** r9 junctions are all axis-aligned, so the stock qiskit-metal junction renderer is correct. The rotated-sheet patch used for the r6–r8 −45° JJs was removed with those designs.

## 2. Coupler geometry comparison (HFSS flux, Q3D Ec)

All runs use the HFSS fixture Q0 + `TC_E` + `TC_N` (SQUIDs 1.41 mm apart), 12 passes, discrete 10–100 MHz, port-only JJ1. Crosstalk = M(neighbour loop ← line) / M(own loop ← own line).

| Variant (`results/flux/…`) | JJ orientation | Own M pH | Crosstalk | Final ΔS | Ec TC (Q3D) |
|---|---|---|---|---|---|
| r7 ring hub, radial SQUID at 45° (`ring_v3_arcT`) | −45° (not Manhattan) | 0.59 / 0.62 | 0.57 % / 1.65 % | 1.9e-4 | 197 MHz |
| ring hub, vertical JJs at a 45° rim notch (`r9_vjj_g44s30`) | 0 / 180° | 0.44 / 0.42 | **16 % / 20 %** | 3.8e-4 | 194 MHz |
| ring hub, radial SQUID facing the flux line (`r9_rad90_o25a50`) | 0 / 90 / 180 / 270° | 0.47 / 0.50 | 2.0 % / 3.5 % | 1.7e-4 | 197 MHz |
| same + necks enter 45 µm above the axis, arc-T 40°/55° (`r9_rad90n_g24s70`) | 0 / 90 / 180 / 270° | 0.62 / 0.63 | 0.74 % / 2.87 % | 2.1e-4 | — |
| rounded-rect hub, vertical-slot SQUID (`rect_w200h150_T`) | 0 / 90 / 180 / 270° | 0.22 / 0.25 | 2.6 % / 2.3 % | 2.9e-4 | 224 MHz |
| **rounded-rect hub, rail SQUID, straight T (`rail_w230h180_T`) → adopted r9** | **0 / 90 / 180 / 270°** | **0.735 / 0.729** | **0.20 % / 0.03 %** | **4.1e-4 (above the 2e-4 target)** | **199 MHz** |

Interpretation:
- **Rail SQUID (best):** a horizontal loop (60 × 20 µm) sits parallel to a straight T, about 30 µm from it. The T's left arm covers the whole loop and the right arm runs under island B, so only one arm couples. The T's symmetric return keeps the far field quadrupole-like.
- **Ring hub:** the necks on the coupler axis limit the arc-T span, so the arms partly cancel. The loop is also about 45 µm from the flux line.
- **Vertical-slot SQUID:** the large islands on both sides screen the narrow window, and only the bottom end of the loop is close to the bar.
- **Vertical JJs in a 45° notch:** the large crescent tab focuses the far field into the loop, which explains the 16–20 % crosstalk.

## 3. Adopted r9 small-model values (Q3D run `…rail_w230h180…`, HFSS run `…rail_w230h180_T…`)

- Qubit: `Ec 193 MHz`, `g_qr 157 MHz`, `2χ 1.38 MHz`, `n_crit 65`, readout↔coupler parasitic `g_cr/g_qr 1.7 %`.
- Coupler: `Ec 199 MHz` (`C_Σ 97.2 fF`), `g_qc 80 MHz`.
- Flux line: own `M 0.73 pH` (`I_Φ0 2.8 mA`), coupler-to-coupler crosstalk `0.20 % / 0.03 %`.

## 4. Limitations

- **Convergence:** the adopted HFSS run ended at ΔS `4.1e-4`, above the `2e-4` target. A 15-pass confirmation is pending.
- **Scope:** all values are single small-model runs (8-pass Q3D; 12-pass HFSS with local return paths only), not full-chip.
- **Not evaluated:** TLS participation of the keyhole gap and of the rounded-rectangle corners.
