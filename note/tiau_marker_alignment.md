# TiAu alignment markers and the wafer die grid

Status (2026-11-19): layout check only. No wafer version with TiAu markers is confirmed for fabrication. The dated
wafer folders, including `20261119/TiAugds/`, are local only (git-ignored). This note records the method and the
tools.

## Constraint

- The TiAu marker photomask `Fluxonium_2Q_bandpass_filter_Mask2.gds` is already fabricated. Its coordinates cannot
  change.
- **Chip markers:** layer 2/0.
  - 4 sites per chip, at (±2529.25, ±4102.5) µm from the chip centre (centre of the 40 µm cross).
  - Each site has a 40 µm cross, two 20 µm crosses and vernier squares.
  - 61 chips, 244 sites.
- **Wafer markers:** layer 3/0.
  - 8 sites. Each has a 2300 × 1800 µm cross, a 550 × 450 µm cross and a label.
  - Big crosses at x = ±22500 µm (y = +36870 and −38110 µm) and at x = ±37350 µm (y = +22620 and −23380 µm).
  - The two halves are not mirror images in y.
- **Chip grid:** the marker mask follows the old single-chip wafer (the laser-cut file).
  - 10 × 10 mm chips on a 10.2 mm pitch (200 µm streets).
  - An odd 9 × 9 grid with a chip centred at (0, 0), rows of 3, 7, 7, 9, 9, 9, 7, 7, 3.

## Finding

The `assemble_wafer.py` grid does not match the marker mask:
- **Pitch:** 10.38 × 10.41 mm (largest die slot plus a 100 µm street).
- **Grid:** an even 8 × 8 grid, with a street crossing at (0, 0).
- **Effect on each die:**
  - Every die is shifted by half a pitch, (5.19, 5.205) mm, from the old chips.
  - The larger pitch then adds 180 µm per column and 210 µm per row.
  - So each chip marker lands somewhere different in every die.
- **20261119 summery:** 37 of the 244 chip-marker sites touch patterns. Some lie across CPW gaps, where TiAu would
  short the centre conductor to ground.

Re-placing the same dies on the 61 old chip centres (`summery_mask2grid`, 61 dies) puts every marker at
(±2507, ±4100) µm in every die (the site-group centre):

| Die | Marker sites |
|---|---|
| D2, D3, D5, D6, Chip4 | All four inside ground (only 5 µm cheese holes nearby) |
| PCM (external Allspice die) | Upper two touch features (a metal pad; a metal edge about 7 µm away); lower two are clear |

**Chip4 and PCM on this grid:**
- They are 10.27–10.28 mm wide, so two of them overlap even when placed diagonally.
- They are placed so that none touches another, and none is on the outermost corner sites.
- The smallest die-to-die gap is 45 µm.
- The outermost die corner is 454 µm inside the wafer edge, the same as the old wafer.

## Tools

- `code/general/tiau_marker_overlay.py`:
  - Copies Mask2 unchanged and writes the marker GDS files (chip, wafer, both, and a wafer view with the 99/0 guide).
  - Overlays them on an assembly (`--stem`, `--assembly-dir`).
  - Classifies each marker site as clear, ground_plane or conflict. A 30 µm window is used, and openings of 5 µm or
    less count as cheese.
- `code/general/assemble_wafer_mask2_grid.py`:
  - Builds the same six die sources as `assemble_wafer.py` on the Mask2 grid.
  - Counts are set in `COUNTS` and the seed is 20261119. Same die-ID style; no street crosses.
  - Checks that no die outlines overlap and that every die fits the wafer.

Not done:
- No dies-and-streets etch artwork exists for the Mask2-grid assembly.
- The PCM features at the two marker sites were not assessed.
- Whether a buried TiAu marker can still be used for alignment depends on the process and was not checked.
