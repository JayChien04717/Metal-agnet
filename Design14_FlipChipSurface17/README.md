# Design14_FlipChipSurface17 — Dual-Tier Flip-Chip Rotated $d=3$ Surface-17

## Stage Status
- **Current stage:** **Stage 1 (Draft layout & figure review)** — **NO Stage 2 EM simulations (Q3D / HFSS / scqubits) or Stage 3 fab GDS have been run yet.**
- **Revision:** `r1` (`2026-11-27`)

## Architecture Summary
`Design14_FlipChipSurface17` splits the `Design13` rotated $d=3$ Surface-17 floating-coupler architecture into two aligned chips bonded face-to-face across a target `~7 µm` vacuum gap via **161 Indium micro-bumps** (Layer 30, `24 × 24 µm` UBM pads on a `350 µm` pitch, `≥ 105 µm` keep-out from all Layer-1 pockets/traces on both tiers):

| Tier | Die Size | Contents | Airbridges |
|---|---|---|---|
| **Q-tier (Top Qubit Chip)** | `8 × 8 mm` | **ONLY** 17 pure circular wide-gap `ScalablePadQubit`s (`256 µm` pad, `46 µm` pocket gap, `make_readout_pad=False`) + 24 `FloatingTunableCoupler`s (`flux_as_ground=True`, no on-chip flux lines) + 161 Indium bump UBM pads | **0** |
| **C-tier (Bottom Interposer Chip)** | `16 × 16 mm` | 56 perimeter launch pads + 4 shared-bus $\lambda/2$ Purcell filters (`6.85–7.18 GHz`) + 17 `VerticalReadoutCoupler`s (`100 µm` circular pads directly beneath qubit islands) feeding 17 exact $\lambda/4$ DRC-B7 readout resonators + 17 vertical XY drive stubs + 24 `VerticalFluxLoop`s directly beneath every coupler SQUID + 161 matching Indium bumps | **217** (C-tier CPW-only; 0 breaking Purcell/readout) |

## Stage 1 DRC Verification (`layout/build.py` exit 0)
- Q-tier control lines / resonators / airbridges: **0 / 0 / 0**
- C-tier pad & `EnclosingCap` collisions: **0**
- C-tier acute (`< 90°`) bends: **0**
- C-tier unbridged Layer-1 CPW crossings: **0** (hard-checked by `_assert_zero_unbridged_cpw_crossings`)
- Short-segment fillet warnings: **0**
- Shared Indium micro-bumps: **161** (`≥ 105 µm` clear of all pockets/traces on both Q-tier and C-tier)

## Entry Points & Outputs
- Build & render script (`qmetal` env): `python -I -X utf8 -B Design14_FlipChipSurface17/layout/build.py`
- Custom C-tier components: `layout/components/flipchip_couplers.py` (`VerticalReadoutCoupler`, `VerticalFluxLoop`)
- Simulation GDS (local):
  - `gds/Design14_2026-11-27_r1_qtier_sim.gds`
  - `gds/Design14_2026-11-27_r1_ctier_sim.gds`
- Stage 1 figures:
  - `fig/qtier.png` — Top Qubit Chip (`8 × 8 mm`)
  - `fig/ctier.png` — Bottom Control & Readout Interposer Chip (`16 × 16 mm`)
  - `fig/chip.png` — Aligned Flip-Chip Stack Overlay
  - `fig/core.png` — 2Q+1TC Flip-Chip Alignment Zoom (`D5 – TC_D5_Z2 – Z2`)

## Unvalidated / Next Stage (Stage 2)
- Vertical flip-chip capacitance ($C_{\text{vert}}$ across `7 µm` gap for readout & XY) and vertical mutual inductance ($M_{fc}$ from `VerticalFluxLoop` to Q-tier coupler SQUID) have not yet been simulated in Q3D / HFSS.
- Parallel-plate chip-mode eigenfrequency check with the 161-bump ground-stitching grid is deferred to Stage 2.
