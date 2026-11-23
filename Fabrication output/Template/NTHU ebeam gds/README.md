# NTHU e-beam: TiAu markers

The markers come from the fixed photomask `Fluxonium_2Q_bandpass_filter_Mask2.gds` (SHA256 `9f0ff079…`, not kept
here). They are stored as two separate references:

| File | Content | Bare rectangle |
|---|---|---|
| `TiAu_marker_wafer_4inch.gds` | 4-inch wafer markers, Mask2 3/0: 8 sites of 2300 × 3155 µm around the wafer, outside the dies | Bounding box + 50 µm (2400 × 3255 µm), no metal |
| `TiAu_marker_chip_single_die.gds` | Chip markers of one die, Mask2 2/0: 4 sites of 84 × 68 µm at (±2507.25, ±4100) µm | Bounding box + 50 µm = 184 × 168 µm ground opening per site (`tiau_ground_opening.py`) |

These files are read in place by the scripts. Do not copy them into a delivery: the delivered
`<date>_<source>_ebeam_NTHU.gds` already carries the markers on `2/0` and `3/0`, and that file plus a README is all
a delivered `NTHU ebeam gds/` folder holds.

Dies must sit on the Mask2 grid: 10 mm dies, 10.2 mm pitch, odd 9 × 9 grid. Rules:
[skill/fabrication_output.md](../../../skill/fabrication_output.md#nthu-tiau-markers).
Background: [note/tiau_marker_alignment.md](../../../note/tiau_marker_alignment.md).
