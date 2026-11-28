# Design13 revision history (Stage 1)

Only `r9` geometry is kept in the code. Older GDS and figures were deleted (2026-11-28); earlier geometries survive in
the run snapshots under `results/*/*/code/` and in the notes listed below. EM simulations were never run on the full chip.

| Revision | Date | Change |
|---|---|---|
| `r4` | 2026-11-27 | 30 mm octilinear Surface-17, circular 256 µm qubits, pill-pad couplers, hand-tuned control routing with oblique (45°) airbridges. |
| `r5` | 2026-11-27 | Clearance-aware octilinear A* router (`code/general/octilinear_router.py`) and checker (`code/general/cpw_clearance.py`). Hard DRC: ≥ 180 µm CPW clearance, 90° crossings only, 35 µm arch, ground-stitching coverage. Rules in `skill/airbridge.md` §1/§4/§4a. |
| `r6` | 2026-11-27 | Circular design language: C-arc couplers, circular SQUID hub, arc flux line. JJs rotated to −45° (later found not Manhattan-compliant). C-arc 36° → 28° after the readout↔coupler Q3D check ([q3d_ro_coupler_crosstalk_2026-11-27.md](q3d_ro_coupler_crosstalk_2026-11-27.md)). |
| `r7` | 2026-11-27 | Qubit Ø280 µm with keyhole readout; ring hub (disc + crescent) with a rim SQUID and arc-T flux line ([coupler_ring_hub_flux_2026-11-27.md](coupler_ring_hub_flux_2026-11-27.md)). |
| `r8` | 2026-11-28 | Trial on a 20 × 20 mm die; meander keep-out hull fix; early sideways turn for crowded readout banks; `u_purcell` 6.10 → 5.60 mm (11/17 resonators on target). |
| `r9` | 2026-11-28 | Back to 30 mm. All JJs Manhattan-compliant (0/90/180/270°). Coupler: rounded-rectangle hub with a rail SQUID and a straight T flux line ([coupler_rect_rail_manhattan_2026-11-28.md](coupler_rect_rail_manhattan_2026-11-28.md)). |
