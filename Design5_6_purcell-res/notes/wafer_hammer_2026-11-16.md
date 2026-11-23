# Design5_6_purcell-res wafer hammer and corner revision

The r3 fabrication candidate keeps the six channel lengths, 400 µm filter-to-resonator
CPW centre spacing, 16/9 µm CPW cross-section and all prime tap positions. The wider
hammer heads shorten the narrow branch on both sides. The larger face gap also
reduces the metal reach into the space between the main CPWs. The user's attached
hammer image supplied the layout reference. Its 50 µm arms and 30 µm head height
describe a different footprint.

| Dimension or measurement | Previous wafer mask | r3 candidate |
|---|---:|---:|
| CPW centre spacing (µm) | 400 | 400 |
| CPW centre to outside of head on each side (µm) | 154 | 125 |
| CPW centre to inside of head on each side (µm) | 179 | 175 |
| Head width along the branch (µm) | 25 | 50 |
| Head height along the main CPW (µm) | 65 | 65 |
| Facing gap (µm) | 42 | 50 |
| Ground clearance at head (µm) | 16 | 16 |
| Local Q3D mutual capacitance (fF) | 1.47064 | 1.47587 |

The mutual capacitance differs by +0.36%. The 6.5 GHz electrostatic comparison uses
two identical 300 µm sections of main CPW, four open ends, the same 500 µm silicon
substrate and common ground. Both cases have eight saved adaptive matrices. The baseline
did not reach the solver convergence threshold and its last two mutual estimates differ
by 0.36%. The candidate reports convergence and its last two mutual estimates differ
by 0.35%. This is a local comparison of mutual capacitance. The full-chip frequency
response and coupling rate have not been revalidated. HFSS was not rerun.

Q3D summary (`results/resonator/20261116_155136_hammer_quick_q3d_ff149ab9/parameter_report.json`)
· Run manifest (`results/resonator/20261116_155136_hammer_quick_q3d_ff149ab9/parameter_report.json`).
The earlier `20261116_154931_hammer_quick_q3d_7b3a6bb7` comparison reused the active
AEDT design and is excluded. The adopted comparison verifies the exact AEDT design
name for each case.

The four cross markers and their windows are removed. The existing lower-right
PFR00 etch is retained. The mask readback checks both launch pads and all six hammer
gaps and heads. The deliverable is positive metal on layer 1/0.

The r3 layout is preserved as an archived fab GDS (`gds/previous/Design5_2026-11-16_wafer_r3_fab.gds`).
The current layout and checks are recorded in the [r5 notes](wafer_hammer_parallel_2026-11-16.md).
