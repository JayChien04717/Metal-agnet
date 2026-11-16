# 10 mm wafer placement candidates — 2026-11-16

All eight layouts use a centred **10.000 × 10.000 mm** die boundary on the qubit layout's positive metal layer **1/0**. The continuous ground perimeter reaches ±5000 µm on both axes. Design2 extends the ground through the former 50 µm edge margin to connect the perimeter to the ground plane.

The current eight-port PCB anchors use x = ±4.5 mm and y = ±2.0 mm on the sides, and x = ±1.1 mm and y = ±4.5 mm on the top and bottom. See [wafer_pcbport.py](../pattern%20design/wafer_pcbport.py). Design4 and Design5 use only the left and right centres, (−4.5, 0) and (4.5, 0) mm. Coordinates are nominal anchors; the launch pad includes the standard inward offset.

| Design | Current fabrication GDS | Geometry companion | Fab mask preview | Perforations filled |
|---|---|---|---|---:|
| Design1 | [Fab](../Design1/gds/Design1_2026-11-16_wafer_r5_fab.gds) | [Unnotched reference](../Design1/gds/Design1_2026-11-16_wafer_r5_reference.gds) | [Preview](../Design1/fig/Design1_2026-11-16_wafer_r5_fab_check.png) | 316 |
| Design2 | [Fab](../Design2/gds/Design2_2026-11-16_wafer_r5_fab.gds) | [Unnotched reference](../Design2/gds/Design2_2026-11-16_wafer_r5_reference.gds) | [Preview](../Design2/fig/Design2_2026-11-16_wafer_r5_fab_check.png) | 207 |
| Design3 | [Fab](../Design3/gds/Design3_2026-11-16_wafer_r5_fab.gds) | [Unnotched reference](../Design3/gds/Design3_2026-11-16_wafer_r5_reference.gds) | [Preview](../Design3/fig/Design3_2026-11-16_wafer_r5_fab_check.png) | 326 |
| Design3_2 | [Fab](../Design3_2/gds/Design3_2_2026-11-16_wafer_r5_fab.gds) | [Unnotched reference](../Design3_2/gds/Design3_2_2026-11-16_wafer_r5_reference.gds) | [Preview](../Design3_2/fig/Design3_2_2026-11-16_wafer_r5_fab_check.png) | 328 |
| Design3_3 | [Fab](../Design3_3/gds/Design3_3_2026-11-16_wafer_r5_fab.gds) | [Simulation geometry](../Design3_3/gds/Design3_3_2026-11-16_wafer_r5_sim.gds) | [Preview](../Design3_3/fig/Design3_3_2026-11-16_wafer_r5_fab_check.png) | 317 |
| Design4 | [Fab](../Design4/gds/Design4_2026-11-16_wafer_r3_fab.gds) | [Simulation geometry](../Design4/gds/Design4_2026-11-16_wafer_r3_sim.gds) | [Preview](../Design4/fig/Design4_2026-11-16_wafer_r3_fab_check.png) | 354 |
| Design5 | [Fab](../Design5/gds/Design5_2026-11-16_wafer_r6_fab.gds) | [Simulation geometry](../Design5/gds/Design5_2026-11-16_wafer_r6_sim.gds) | [Preview](../Design5/fig/Design5_2026-11-16_wafer_r6_fab_check.png) | 32 |
| Design6 | [Fab](../Design6/gds/Design6_2026-11-16_wafer_r7_fab.gds) | [Simulation geometry](../Design6/gds/Design6_2026-11-16_wafer_r7_sim.gds) | [Preview](../Design6/fig/Design6_2026-11-16_wafer_r7_fab_check.png) | 315 |

The exporter reads each output back to check the continuous metal perimeter and the signal and slots at the pad, taper, and neck of every launch. Design2 uses straight lead-ins from the standard pad coordinates. Each active port uses one rounded turn on the new section to match the direction of the retained route. Both Design2 and Design6 ground the spare P3 at its neck.

Design2 engraves [FTFC00](../Design2/fig/Design2_2026-11-16_wafer_die_id.png) at the lower right. Design4 and Design5 omit corner test junctions. The current Design2 and Design3_3 fab masks use F at the upper right and lower left, and T at the upper left and lower right. Design6 uses T at all four corners.

The masks place the vertical coupons at x = ±4.20 mm and y = ±4.00 mm. On Design2, Design3_3 and Design6 the wafer F/T marks measure 260 × 350 µm with 55 µm strokes. Design3_3 and Design6 replace the old cross polygons and the surrounding 5 µm margin with one ground polygon at each corner. Design2's parent mask had no crosses. The die ID etches are retained. The fluxonium coupon's 60 µm array envelope clears ground by 55 µm.

The previous fab masks were reproduced in Qiskit Metal and matched by KLayout XOR before the keep-out repair. The repair fills perforations that extend into the outer 500 µm or sit within 50 µm of a device etch. The fabrication XOR contains only those added ground areas. The selected CPW routes, device contacts, test structures, and layout parameters follow the preceding design revision. The simulation and unnotched companions were exported from that geometry with their documented layer conventions.

Use the qmetal environment and KLayout. Run `python -B code/general/export_wafer_test_junctions.py --design DesignN --output reproduced_fab.gds` to reproduce a fab mask and add `--reference-only` with a separate output filename for its geometry companion. The coupon pair uses [the DesignTestJunction entry points](../DesignTestJunction/README.md).

The four local fabrication handoff files in `20261116` match the Design copies by SHA-256. The local assembly `20261116/summery.gds` uses those masks. The dated handoff folders are outside this repository update. The parent fab masks are in each Design's `gds/previous/` folder. No HFSS or Q3D solve was run for the fabrication keep-out repair.
