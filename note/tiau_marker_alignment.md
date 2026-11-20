# NTHU TiAu alignment markers and the 4-inch wafer grid

Self-contained geometry and alignment reference for the pre-fabricated **NTHU TiAu mask** (e-beam alignment-marker photomask), written so any reader can reproduce or check wafer alignment **without** needing the original photomask GDS file.

Operational delivery rules, SHA256 hashes, and pre-tapeout checklist live in [skill/fabrication_output.md](../skill/fabrication_output.md).

---

## 1. What the NTHU TiAu mask is

In the original 4-inch (100 mm) wafer fabrication set, two companion GDS files shared the exact same 9 × 9 wafer layout:
1. **Original base-metal wafer** (`Fluxonium_2Q_bandpass_filter_送雷切.gds`): 61 die slots (`10 × 10 mm` each) on a `10.2 mm` pitch (`200 µm` dicing streets, spanning `±45.9 mm`), with base metal on layers `1/10`, `1/11`, `1/100` plus alignment markers on `2/0` and `3/0`.
2. **NTHU TiAu mask** (original file `Fluxonium_2Q_bandpass_filter_Mask2.gds`): the lift-off photomask derived from that same wafer by dropping Layer 1 and keeping **only** the TiAu alignment markers:
   - **Layer `2/0` (chip markers):** 2,196 polygons = 61 die slots × 4 sites/die × 9 polygons/site.
   - **Layer `3/0` (wafer markers):** 24 polygons = 8 wafer-level alignment sites.

Because this physical TiAu photomask is already fabricated and reused for NTHU e-beam runs, its marker coordinates are fixed. You do **not** need the original photomask file to work with this repository: the exact marker polygons are extracted unchanged into two tracked reference files under `Fabrication output/Template/NTHU ebeam gds/`:
- `TiAu_marker_wafer_4inch.gds` (layer `3/0`, 8 wafer sites)
- `TiAu_marker_chip_single_die.gds` (layer `2/0`, 4 chip sites of one die centred at `(0, 0)`)

---

## 2. Exact wafer grid, populated slots, and marker coordinates

Every grid slot sits at `(x, y) = (ix × 10.2 mm, iy × 10.2 mm)` on an **odd 9 × 9 grid** (`ix, iy ∈ {−4, …, +4}`) with a die centred at `(0, 0)`.

* **Photomask chip-marker slots (61 total):** row counts `3, 7, 7, 9, 9, 9, 7, 7, 3` (`X` + `o` below).
* **Populated dies for wafer assembly (49 dies, inner 7 × 7):** **omit the outermost 3-chip strip on each of the four sides** (`|ix| = 4` and `|iy| = 4`, marked `o` below, 12 slots total), because wafer clamping presses on those edge chips during processing. Populate only the inner **7 × 7 = 49 slots** (`ix, iy ∈ {−3, …, +3}`, marked `X`, spanning `±35.7 mm`):

```text
iy = +4 (y = +40.8 mm):   .  .  .  o  o  o  .  .  .   <- omit (wafer clamp)
iy = +3 (y = +30.6 mm):   .  X  X  X  X  X  X  X  .   (7 populated)
iy = +2 (y = +20.4 mm):   .  X  X  X  X  X  X  X  .   (7 populated)
iy = +1 (y = +10.2 mm):   o  X  X  X  X  X  X  X  o   (7 populated; omit ix = ±4)
iy =  0 (y =   0.0 mm):   o  X  X  X  X  X  X  X  o   (7 populated; omit ix = ±4)
iy = -1 (y = -10.2 mm):   o  X  X  X  X  X  X  X  o   (7 populated; omit ix = ±4)
iy = -2 (y = -20.4 mm):   .  X  X  X  X  X  X  X  .   (7 populated)
iy = -3 (y = -30.6 mm):   .  X  X  X  X  X  X  X  .   (7 populated)
iy = -4 (y = -40.8 mm):   .  .  .  o  o  o  .  .  .   <- omit (wafer clamp)
ix =                     -4 -3 -2 -1  0 +1 +2 +3 +4
```

* **Die size & streets:** each populated `X` is a `10.0 × 10.0 mm` die (`±5.0 mm` local), leaving `200 µm` dicing streets between adjacent dies (`10.2 mm` pitch).
* **Chip markers (layer `2/0`, inside every die):** 4 sites per die, each an `84 × 68 µm` group (one `40 µm` cross, two `20 µm` crosses, and vernier squares) centred at `(±2507.25, ±4100) µm` from the die centre (40 µm cross centre at `(±2529.25, ±4102.5) µm`). Each site requires a `184 × 168 µm` bare-substrate opening in base metal (`+50 µm` margin per side).
* **Wafer markers (layer `3/0`, outside the dies):** 8 sites sitting in the unoccupied corner regions of the 9 × 9 grid, each a `2300 × 3155 µm` group (one `2300 × 1800 µm` cross, one `550 × 450 µm` cross, and a label) centred at `(±22500, +37547.5)`, `(±22500, −37432.5)`, `(±37350, +23297.5)`, and `(±37350, −22702.5) µm` (not mirror-symmetric in y). Each requires a `2400 × 3255 µm` bare rectangle (`+50 µm` margin per side, clearing the nearest street edge by 55 µm).

---

## 3. Why the generic 8 × 8 grid conflicts with the NTHU TiAu mask

The generic `code/general/assemble_wafer.py` grid uses a `10.38 × 10.41 mm` pitch on an **even 8 × 8 grid** (street crossing at `(0, 0)`, no die at the origin):
- **Half-pitch origin shift:** every die is offset by `(5.19, 5.205) mm` relative to the NTHU TiAu grid.
- **Cumulative pitch drift:** `180 µm/column` in x and `210 µm/row` in y.
- **Result:** chip markers land at a different local position on every die and can cross CPW gaps, shorting centre conductors to ground with TiAu metal.

Placing dies on the inner 7 × 7 centres of Section 2 locks all four chip-marker sites to `(±2507.25, ±4100) µm` in every die, where standard 10 mm dies have solid ground plane (only 5 µm cheese holes nearby). Slightly oversized external dies (`~10.27–10.28 mm`, wider than the `10.20 mm` pitch) must be placed so none are adjacent, preserving die-to-die clearance (`≥ 45 µm`).

---

## 4. Shared alignment tools (`code/general/`)

- **`tiau_marker_overlay.py`:** extracts chip/wafer markers, overlays them on an assembly (`--stem`, `--assembly-dir`), and classifies every site as `clear`, `ground_plane`, or `conflict` (30 µm check window; ≤ 5 µm openings treated as cheese).
- **`assemble_wafer_mask2_grid.py`:** places dies on the NTHU TiAu mask grid (note: verify the script omits the 12 outermost `|ix| = 4` / `|iy| = 4` clamp slots before generating a new NTHU wafer delivery), checking that no die outlines overlap and every die fits inside the 4-inch wafer.
- **`tiau_ground_opening.py`:** cuts `184 × 168 µm` bare-substrate rectangles (`--margin 50`) in layer `1/0` at all four chip-marker sites per flat source cell of an unrotated NTHU TiAu grid assembly, verifies all sites are clear and no die ID is touched, and records the openings in `<out-stem>.json`.

---

## 5. Open items

- Dies-and-streets etch artwork for the NTHU TiAu grid assembly is not yet automated.
- Update `assemble_wafer_mask2_grid.py` if needed so it populates only the inner 7 × 7 (49 dies) and omits the 12 outermost clamp slots (`|ix| = 4` or `|iy| = 4`).
- Whether a buried TiAu marker (covered by base metal without an opening) can still be used for alignment is process-dependent and unverified; bare-substrate rectangles are therefore required.
