# Design9_CPW_Bandpass_Filter figures

Stage 1 draft layout figures rendered directly from Qiskit Metal QGeometry (with exact circular-arc
fillets) by [`layout/build.py`](../layout/build.py). **No GDS file has been generated yet**
(per user instruction; awaiting layout review).

- Grey is metal and white is bare silicon substrate.
- Revision: `2026-11-21_r1`

| Figure | Shows | Window x0, x1, y0, y1 (µm) |
|---|---|---|
| [chip.png](chip.png) | Whole chip (10 × 6 mm) | -5000, 5000, -3000, 3000 |
| [launch_pad.png](launch_pad.png) | LP_IN + feed_in + T1 | -4980, -2800, -650, 650 |
| [t_junction_stub.png](t_junction_stub.png) | T3 + folded λ/4 stub S3 | 281.25, 2418.75, -250, 2600 |
| [connecting_loop.png](connecting_loop.png) | λ/4 connecting U-loop L23 (T2–T3) | -1250, 1250, -1800, 320 |
| [tap_closeup.png](tap_closeup.png) | T1 & S1_short (10/6 µm CPW) | -3180, -2200, -180, 520 |
