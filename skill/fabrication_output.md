# Fabrication output: 4-inch wafer assembly, streets, die IDs, tones, and e-beam markers

Every GDS that leaves the repository for fabrication is written under `Fabrication output/`, grouped first by **date** (`<YYYYMMDD>/`), then by **tool / e-beam site**.

Single-die mask export (fab swap, launch pads, cheesing, test coupons) is in [gds_export.md](gds_export.md). E-beam junction mask sizing and drawing is in [ebeam_junction.md](ebeam_junction.md).

---

## 1. Required content of every delivered wafer GDS

Every laser-writer wafer delivery must contain all four items below:

| # | Item | Rule | Tool |
|---|---|---|---|
| 1 | **Dicing streets** | Street lines on the delivery's die grid, extended to the wafer boundary, identical in positive and etch files (NTHU TiAu grid: 10 mm dies, 10.2 mm pitch, 200 µm streets; standard 8 × 8 grid: 100 µm streets) | `assemble_wafer.py` / `assemble_wafer_mask2_grid.py`, `invert_wafer_tone.py` |
| 2 | **Die IDs** | Unique consecutive ID per populated die (`D00`, `D01`, …) in the upper-right corner, distinct from the source design ID (Section 1.1) | Assembler |
| 3 | **Positive and negative (etch)** | Two files from **one** positive source: `*_positive.gds` (metal = polygons) and `*_etch.gds` (exposed = polygons), plus `*_etch.dxf` if required. Etch must be the exact complement inside the exposed scope; verify text polarity after inversion | `invert_wafer_tone.py --source-stem …` |
| 4 | **E-beam site alignment markers** | Site-required markers in **bare rectangles** (no metal within marker bbox + clearance margin). Modify the positive source **before** inverting so both tones agree | Section 2 |

Junction layers never go in laser-writer files; they go in the site's e-beam file (Section 3).

### 1.1 4-inch wafer placement IDs and tone inversion

<a id="4-inch-wafer-placement-ids"></a>

* **Ordering:** sort occupied die positions by descending y, then ascending x (left to right, top to bottom), skipping empty wafer-edge slots. Number consecutively from `D00` (≥ 2 digits, e.g. `D00`–`D51` for a 52-die wafer).
* **Glyphs:** polygon text on the process layer, 100 µm nominal height, right/top bounds 50 µm inside the source die edges. Etch into solid ground margin (or positive metal where the vendor-die corner is blank). Check the label window for collisions before writing; preserve street width and alignment crosses.
* **Assembly record:** record each die's ID, cell, region, centre, and local glyph bounds in the assembly JSON. Use shared source cells + one wrapper per die to keep repeated electrical geometry identical.
* **Tools:** `code/general/assemble_wafer.py` (with optional `--output-stem`) builds the positive-metal GDS/JSON/PNG; `code/general/invert_wafer_tone.py --source-stem …` generates the matching etch GDS (`1/0` exposure, `99/0` wafer guide) and closed-polyline µm DXF (`L1D0` / `L99D0`). Read both back and verify the `D00`–last sequence, upper-right positions, street geometry, and text polarity.

Exposure scopes of `invert_wafer_tone.py`:

| Option | Exposed | Kept metal | Output suffix |
|---|---|---|---|
| `--dies-and-streets` (default) | Die etch + street lines extended to wafer boundary | Alignment crosses, rest of wafer | `_etch_dies_and_streets` (+ verification JSON) |
| `--dies-only` | Die etch only | Everything outside die rectangles | `_etch_dies_only` |
| `--dies-and-crosses` | Die etch + street alignment crosses (250 × 250 µm, 5 µm strokes) | Streets, wafer perimeter | `_etch_dies_and_crosses` |
| `--wafer` | Everything inside wafer guide outer edge (alignment crosses stay retained metal) | — | `_etch` |

---

## 2. Alignment markers by e-beam site

| Site | Markers | Status |
|---|---|---|
| **NTHU** | TiAu markers: 4-inch wafer markers + chip markers (below) | Defined |
| **TSRI** | TSRI A12 marks (`A12mark_TSRI.gds`) | **Pending** (placement, polarity, layer TBD; no TSRI delivery until defined) |
| **NCKU** | None | Base unchanged |

### NTHU: TiAu markers

NTHU e-beam runs reuse an already-fabricated **NTHU TiAu mask** (lift-off photomask carrying only layers `2/0` chip markers and `3/0` wafer markers for an odd 9 × 9, 10.2 mm-pitch wafer grid). Because its coordinates are fixed, the base metal must fit around them. **Omit the outermost 3-chip strip on each of the four sides** (`|ix| = 4` and `|iy| = 4`, 12 slots total) because wafer clamping presses on them, leaving the inner **7 × 7 = 49 populated dies**. Readers do not need the original photomask file: the exact marker polygons are stored in `Fabrication output/Template/NTHU ebeam gds/`, and the full grid map and coordinates are documented in [note/tiau_marker_alignment.md](../note/tiau_marker_alignment.md):

| | **4-inch wafer markers** | **Chip markers** (inside each die) |
|---|---|---|
| Reference (`Template/NTHU ebeam gds/`) | `TiAu_marker_wafer_4inch.gds` (NTHU TiAu mask 3/0) | `TiAu_marker_chip_single_die.gds` (NTHU TiAu mask 2/0, die at (0, 0)) |
| Sites | 8 sites: 2300 × 1800 µm cross, 550 × 450 µm cross, label | 4 sites/die: 40 µm cross, two 20 µm crosses, vernier squares |
| Site bounding box | 2300 × 3155 µm at `(±22500, +37547.5)`, `(±22500, −37432.5)`, `(±37350, +23297.5)`, `(±37350, −22702.5)` µm | 84 × 68 µm at `(±2507.25, ±4100)` µm from die centre (40 µm cross at `(±2529.25, ±4102.5)`) |
| Bare rectangle | Bbox + 50 µm/side (2400 × 3255 µm); no die, die ID, street, or metal inside (closest street 55 µm away) | Bbox + 50 µm/side = 184 × 168 µm opening in ground metal (`tiau_ground_opening.py --margin 50`) |
| Grid & compliance | NTHU TiAu mask grid (10 mm dies, 10.2 mm pitch, inner 7 × 7 = 49 populated dies with die at (0, 0); omit outermost 3-chip clamp strips `|ix|=4` / `|iy|=4`; `assemble_wafer_mask2_grid.py`) | Same grid; `tiau_ground_opening.py` cuts all 4 openings per source cell |
| Check (`tiau_marker_overlay.py`) | 8 / 8 `outside_dies/clear` | 4 / die `clear`, no metal in opening, no die ID touched |

* **Bare rectangle polarity:** empty in positive, fully exposed in etch (including wafer-marker rectangles outside dies).
* Any wafer not on the NTHU TiAu mask grid cannot go to NTHU.

### TSRI: A12 marks (pending)

Reference: `Template/TSRI ebeam gds/A12mark_TSRI.gds` (12 × 12 mm field on layer 0/0: four 4 × 4 fine-cross groups at `{±3100…±3400} µm`, four 1200 µm coarse crosses at `(±3600, ±4200) µm`, two 500 × 50 µm diagonal corner L marks at `(±6000, ±6000) µm`). Before any TSRI delivery, define polarity, target layer, per-die vs. per-field placement, bare-rectangle margin, and wafer-level mark set.

---

## 3. E-beam GDS (per site)

- **Filename:** `<SITE> ebeam gds/<date>_<source>_ebeam_<SITE>.gds`.
- **Junction layers:** 1–4, or 11–14 when base sits on 1/0 ([ebeam_junction.md](ebeam_junction.md)); state the layer map in README and manifest.
- **Base metal:** reference layer only, vertex-identical to the laser-writer positive of the same delivery (matching SHA256).
- **Markers:** include site marker artwork on its reference layers (NTHU: chip `2/0`, wafer `3/0`).
- **Sizing:** from `code/qubit/jj_resistance_model.py`; record EJ inputs and calibration in README.

---

## 4. Folder layout and marker references

```text
Fabrication output/
  Template/                          tracked: empty site folders + marker references + README
    Laser writer gds/
    NTHU ebeam gds/   TiAu_marker_wafer_4inch.gds               4-inch wafer markers (3/0)
                      TiAu_marker_chip_single_die.gds           chip markers of one die (2/0)
    TSRI ebeam gds/   A12mark_TSRI.gds                          pending
    NCKU ebeam gds/   (no marker)
  <YYYYMMDD>/                        one delivery; local only (git-ignored)
    README.md, manifest.json
    Laser writer gds/                <date>_<source>_positive.gds, _etch.gds (+ _etch.dxf)
    <SITE> ebeam gds/                <date>_<source>_ebeam_<SITE>.gds
```

- Copy `Template/` to `<YYYYMMDD>/` and delete unused site folders. Never edit `Template/` or overwrite a delivered folder (revisions use `<YYYYMMDD>_r2` and name what they replace).
- `<source>`: wafer stem or `<Design>_<rev>` (no spaces).

| Fixed marker reference (never edit) | SHA256 |
|---|---|
| `NTHU ebeam gds/TiAu_marker_wafer_4inch.gds` (NTHU TiAu mask 3/0, 24 polys) | `7eafe6bec85db93ff614aa0e705673bcb0cdebd43cc3527b3315acddfe989f78` |
| `NTHU ebeam gds/TiAu_marker_chip_single_die.gds` (NTHU TiAu mask 2/0 ±5 mm, 36 polys) | `d66472dbaec88d81aa89b861ee7e23d2ec0519a2b73eac36d684c76958ab98fb` |
| `TSRI ebeam gds/A12mark_TSRI.gds` | `8d91a98d5b449c86eff7c0abfe6f32dea42ffd85907c6e938ab69a2629149aad` |

---

## 5. Delivery README and `manifest.json`

- **`README.md` (English, concise):** date, revision, replaced folder, e-beam site, source design/revision + `final_selection` status, tone & layer map per file, street grid & die-ID range, markers & margin, checks done, and pending/unvalidated items.
- **`manifest.json`:** record `date`, `revision`, `replaces`, `ebeam_site`, `source` (`gds`, `sha256`, `design_revision`, `final_selection`), `grid`, `markers`, `files` (each with `path`, `sha256`, `tone`/`scope`/`made_from`, `layers`, `ej_inputs_ghz`), `scripts`, `checks`, and `pending`.

---

## 6. Checklist before sending

1. **Streets:** present, extended to wafer edge, identical in positive and etch.
2. **Die IDs:** unique, consecutive from `D00`, upper-right corner, verified polarity in both tones.
3. **Tones:** positive and etch from one source, complement check passed, DXF matches GDS.
4. **Markers:** NTHU on NTHU TiAu mask grid (inner 7 × 7 = 49 populated dies; outermost 3-chip clamp strips omitted) with 8/8 wafer + all chip sites clear and bare rectangles empty in positive / exposed in etch; TSRI blocked while pending; NCKU none.
5. **E-beam file:** holds junction layers + marker references; base reference matches laser-writer positive.
6. **Units & provenance:** µm (1 nm precision), new dated folder, all SHA256 hashes in `manifest.json`, unconfirmed items listed in `pending`.
7. **Confirmation:** source version explicitly user-confirmed ([parameter_provenance.md §6](parameter_provenance.md#6-final-version-confirmation-before-commitpush)).

*Current tooling gap:* no single-step delivery script yet, and the NTHU TiAu mask grid wafer still lacks dies-and-streets etch artwork ([note/tiau_marker_alignment.md](../note/tiau_marker_alignment.md)).
