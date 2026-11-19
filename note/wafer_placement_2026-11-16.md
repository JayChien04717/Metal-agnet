# 10 mm wafer placement candidates — 2026-11-16

## Transmon flux-ground candidates

The 18 November candidates meet the 4 µm minimum between the SQUID pocket and the adjacent flux-line slot. The values below are shortest distances measured on the written positive-metal GDS at 1 nm database resolution. The fabrication and companion variants give the same minimum at each device.

| Design and device | Parent minimum (µm) | Candidate minimum (µm) | Candidate fabrication mask |
|---|---:|---:|---|
| Design2 TC | 1.591 | 5.091 | r9 (`Design2_2026-11-18_wafer_r9_fab.gds`, removed 2026-11-18; in git history) |
| Design3 TC | 3.000 | 4.500 | r8 (`Design3_3_2026-11-18_wafer_r8_fab.gds`, removed 2026-11-18; in git history) |
| Design6 Q1, Q2, TC | 3.000 each | 4.500 each | r11 TXT01 (`Design6_2026-11-18_wafer_r11_txt01_fab.gds`, removed 2026-11-18; in git history) |

The existing Design1, former Design3 and former Design3_2 TC masks and their unnotched companions each measure 4.000 µm. The candidates await final version confirmation. No HFSS or Q3D solver was rerun.

From the repository root, select Design2, Design3 or Design6 and use unused output paths. The exporter reads the geometry inputs and regional bounds from the selected design.json and checks the width after writing each GDS.

```powershell
conda run -n qmetal python -B code/general/export_wafer_test_junctions.py --design Design6 --flux-ground-review --output Design6_flux_review_fab.gds
conda run -n qmetal python -B code/general/export_wafer_test_junctions.py --design Design6 --flux-ground-review --reference-only --output Design6_flux_review_sim.gds
```

## Fluxonium charge-ground candidate

Design3 r9 builds on the r8 TC flux-ground candidate. Retracting the FX2 charge-line tip by 1 µm raises its minimum ground strip from 2.000 to 3.000 µm. FX1 measures 3.000 µm on both revisions. The candidate awaits final confirmation. No HFSS or Q3D solver was rerun. Fabrication GDS (`Design3_3_2026-11-18_wafer_r9_fab.gds`, removed 2026-11-18; in git history) · Simulation GDS (`Design3_3_2026-11-18_wafer_r9_sim.gds`, removed 2026-11-18; in git history)

The fabrication masks for Design1, former Design3 and former Design3_2 give 9.000 µm at both fluxonium charge-line tips. The Design2 candidate gives 12.500 µm at both tips.

```powershell
conda run -n qmetal python -B code/general/export_wafer_test_junctions.py --design Design3 --charge-ground-review --output Design3_charge_review_fab.gds
conda run -n qmetal python -B code/general/export_wafer_test_junctions.py --design Design3 --charge-ground-review --reference-only --output Design3_charge_review_sim.gds
```

## Existing wafer layouts

All eight layouts use a centred **10.000 × 10.000 mm** die boundary on the qubit layout's positive metal layer **1/0**. The continuous ground perimeter reaches ±5000 µm on both axes. Design2 extends the ground through the former 50 µm edge margin to connect the perimeter to the ground plane.

The current eight-port PCB anchors use x = ±4.5 mm and y = ±2.0 mm on the sides, and x = ±1.1 mm and y = ±4.5 mm on the top and bottom. See [wafer_pcbport.py](../pattern%20design/wafer_pcbport.py). Design4 and Design5 use only the left and right centres, (−4.5, 0) and (4.5, 0) mm. Coordinates are nominal anchors; the launch pad includes the standard inward offset.

| Design | Current fabrication GDS | Geometry companion | Fab mask preview | Ground openings |
|---|---|---|---|---:|
| Design1 | [Fab](../Design1/gds/Design1_2026-11-16_wafer_r7_fab.gds) | [Unnotched reference](../Design1/gds/Design1_2026-11-16_wafer_r7_reference.gds) | Preview (`Design1_2026-11-16_wafer_r7_fab_check.png`, removed 2026-11-19; in git history) | 174,864 |
| Design2 | [Fab](../Design2/gds/Design2_2026-11-16_wafer_r8_fab.gds) | [Unnotched reference](../Design2/gds/Design2_2026-11-16_wafer_r8_reference.gds) | Preview (`Design2_2026-11-16_wafer_r8_fab_check.png`, removed 2026-11-19; in git history) | 171,388 |
| former Design3 | removed 2026-11-18 (git history: `Design3/gds/Design3_2026-11-16_wafer_r7_fab.gds`) | | | 174,597 |
| former Design3_2 | removed 2026-11-18 (git history: `Design3_2/gds/Design3_2_2026-11-16_wafer_r7_fab.gds`) | | | 174,603 |
| Design3 (former Design3_3) | [Fab](../Design3/gds/Design3_3_2026-11-16_wafer_r7_fab.gds) | [Simulation geometry](../Design3/gds/Design3_3_2026-11-16_wafer_r7_sim.gds) | Preview (`Design3_3_2026-11-16_wafer_r7_fab_check.png`, removed 2026-11-19; in git history) | 173,429 |
| Design4 | [Fab](../Design4/gds/Design4_2026-11-16_wafer_r4_fab.gds) | [Simulation geometry](../Design4/gds/Design4_2026-11-16_wafer_r4_sim.gds) | Preview (`Design4_2026-11-16_wafer_r4_fab_check.png`, removed 2026-11-19; in git history) | 192,246 |
| Design5 | [Fab](../Design5/gds/Design5_2026-11-16_wafer_r7_fab.gds) | [Simulation geometry](../Design5/gds/Design5_2026-11-16_wafer_r7_sim.gds) | Preview (`Design5_2026-11-16_wafer_r7_fab_check.png`, removed 2026-11-19; in git history) | 175,954 |
| Design6 | [Fab r10 TXT01](../Design6/gds/Design6_2026-11-17_wafer_r10_txt01_fab.gds) | [Simulation geometry r10](../Design6/gds/Design6_2026-11-17_wafer_r10_sim.gds) | Preview (`Design6_2026-11-17_wafer_r10_txt01_fab_check.png`, removed 2026-11-19; in git history) | 160,422 |

The exporter reads each output back to check the continuous metal perimeter and the signal and slots at the pad, taper, and neck of every launch. Design2 uses straight lead-ins from the standard pad coordinates. Each active port uses one rounded turn on the new section to match the direction of the retained route. Both Design2 and Design6 ground the spare P3 at its neck.

The six designs that carry etched die IDs place them in the upper-left ground margin, 50 µm inside both edges. Design2 uses `FTFC00`. Design4 and Design5 omit corner test junctions. The current Design2 and Design3 fab masks use F at the upper right and lower left, and T at the upper left and lower right. Design6 uses T at all four corners.

The masks place the vertical coupons at x = ±4.20 mm and y = ±4.00 mm. On Design2, Design3 and Design6 the wafer F/T marks measure 260 × 350 µm with 55 µm strokes. Design3 and Design6 replace the old cross polygons and the surrounding 5 µm margin with one ground polygon at each corner. Design2's parent mask had no crosses. The die ID etches are retained. The fluxonium coupon's 60 µm array envelope clears ground by 55 µm.

The current Fab GDS uses 5 × 5 µm holes at a 20 µm centre pitch. KLayout measured all four outermost hole edges at ±4,452.5 µm for each 10 mm die. This matches the Design3 reference margin of 547.5 µm. Per-hole clearance from other etched features is at least 50 µm. The geometry companions preserve their original solver or unnotched geometry.

Use the qmetal environment and KLayout. Run `python -B code/general/export_wafer_test_junctions.py --design DesignN --output reproduced_fab.gds` to reproduce a fab mask and add `--reference-only` with a separate output filename for its geometry companion. The coupon pair uses [the DesignTestJunction entry points](../DesignTestJunction/README.md).

The current grids follow the [Design3 launch profile and four-edge rule](../skill/gds_export.md#3-checks-after-an-export). The hole sizes, pitch and bounds were checked by reading the written GDS. The side-port comparison (`Design3_3_2026-11-16_wafer_r7_cheese_edge_check.png`, removed 2026-11-19; in git history) shows the same boundary in Design2 and Design3. The earlier 50 µm-grid Design2 launch restoration is documented in its r6 to r7 preview (`Design2_2026-11-16_wafer_r7_launch_cheese_check.png`, removed 2026-11-19; in git history). No solver was rerun for the pitch change.

The 17 November assembly (local folder `20261117`, deleted on 2026-11-19, never tracked) carried the Design6 r10 electrical geometry and the `TXT01` fabrication ID. The [16 November assembly](../20261116/summery.gds) contains the Design6 r9 parent (`Design6_2026-11-16_wafer_r9_fab.gds`, removed 2026-11-18; in git history). The original [r10 tail update](../Design6/README.md#loaded-readout-s21) passed the KLayout comparison outside its four tail windows. No HFSS or Q3D solve was run for the fabrication keep-out repair, the r10 GDS export or the die-ID update.

The companions `20261116/summery_etch.gds` and `20261116/summery_etch.dxf` describe the exposure across the full wafer inside the outer edge of the wafer guide. The launch slots, cheese holes, open streets and blank wafer perimeter all expose. The alignment crosses are retained-metal openings in the exposure. In GDS use `1/0` for exposure and `99/0` for the guide. The DXF uses closed polylines in micrometres, with `L1D0` for exposure and `L99D0` for the guide. Run `python -B code/general/invert_wafer_tone.py 20261116 --wafer` with KLayout available to reproduce this historical whole-wafer variant. The default export now uses 100 µm streets and retained metal crosses.

The `summery_etch_dies_only.gds` and `summery_etch_dies_only.dxf` variants confine exposure to the 52 die rectangles. They carry a single geometry layer (`1/0` in GDS and `L1D0` in DXF). Use µm when importing the DXF. Run the same command with `--dies-only` to reproduce the pair.

The `summery_etch_dies_and_crosses.gds` and `summery_etch_dies_and_crosses.dxf` variants also expose the 57 original street alignment crosses. Each etched cross spans 250 × 250 µm with 5 µm strokes. The rest of the streets and the wafer perimeter retain metal. Run the same command with `--dies-and-crosses` to reproduce the pair.
