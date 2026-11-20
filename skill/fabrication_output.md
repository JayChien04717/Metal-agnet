# Fabrication output: delivered GDS, streets, die IDs, tones and e-beam markers

Every GDS that leaves the repository for fabrication is written under `Fabrication output/`. A delivery is grouped first
by **date**, then by **tool / e-beam site**.

How the geometry itself is exported is covered elsewhere:
- [gds_export.md §2.1](gds_export.md#21-fab-gds-workflow): the fab GDS workflow. The simulation parameters are
  confirmed first, then the design is swapped to the fab components, the user is asked about the site, markers and
  junctions, and cheese is added last.
- [gds_export.md](gds_export.md): masks, cheesing, die-ID style and street geometry.
- [ebeam_junction.md](ebeam_junction.md): junction layers. This skill defines what
every delivered file must contain and how it is packaged.

## 1. Required content of every delivered wafer GDS

Every laser-writer file must have all four items below. A file missing any of them is not a deliverable.

| # | Item | Rule | Tool |
|---|---|---|---|
| 1 | **Dicing streets** | Street lines on the delivery's die grid, extended to the wafer boundary. They are identical in the positive and etch files (on the Mask2 grid: 10 mm dies, 10.2 mm pitch, 200 µm streets) | `assemble_wafer.py` / `assemble_wafer_mask2_grid.py`, `invert_wafer_tone.py` |
| 2 | **Die IDs** | One unique ID per populated die (D00, D01, …), in the upper-right corner, ordered as in [gds_export.md](gds_export.md#4-inch-wafer-placement-ids). Kept separate from the design identifier | Assembler |
| 3 | **Positive and negative (etch)** | Two files from **one** positive source: `*_positive.gds` (metal = polygons) and `*_etch.gds` (exposed = polygons). Add `*_etch.dxf` if the writer needs it. The etch file must be the exact complement of the positive inside the exposed scope, and the text polarity must be checked after inversion | `invert_wafer_tone.py --source-stem …` |
| 4 | **Alignment markers for the e-beam site** | The markers required by the site that will write the junctions (Section 2). Each marker sits in a **bare rectangle**: no metal within the marker bounding box grown by the clearance margin. The base metal is changed in the positive source **before** inverting, so both tones agree | Section 2 |

Junction layers are never in laser-writer files. They go in the e-beam file of the site (Section 3).

## 2. Alignment markers by e-beam site

| Site | Markers | Status |
|---|---|---|
| **NTHU** | TiAu markers, two classes: 4-inch wafer markers and chip markers (below) | Defined |
| **TSRI** | TSRI A12 marks (`A12mark_TSRI.gds`) | **Pending**: placement, polarity and layer to be defined with TSRI later. No TSRI-bound delivery until then |
| **NCKU** | None | Base unchanged |

### NTHU: TiAu markers

The TiAu markers come from an already-fabricated photomask, `Fluxonium_2Q_bandpass_filter_Mask2.gds`, so their
coordinates cannot change. The base metal must fit around them:

| | **4-inch wafer markers** (around the wafer) | **Chip markers** (inside each die) |
|---|---|---|
| Reference file (`Template/NTHU ebeam gds/`) | `TiAu_marker_wafer_4inch.gds` (Mask2 layer 3/0) | `TiAu_marker_chip_single_die.gds` (Mask2 layer 2/0, die at (0, 0)) |
| Sites | 8 sites, each a 2300 × 1800 µm cross, a 550 × 450 µm cross and a label | 4 sites per die, each a 40 µm cross, two 20 µm crosses and vernier squares |
| Site bounding box | 2300 × 3155 µm, centred at (±22500, +37547.5), (±22500, −37432.5), (±37350, +23297.5) and (±37350, −22702.5) µm | 84 × 68 µm, centred at (±2507.25, ±4100) µm from the die centre (40 µm cross at (±2529.25, ±4102.5)) |
| Where it sits | Outside the dies, between die rows near the wafer edge | Inside every die |
| Bare rectangle | Bounding box + 50 µm per side (2400 × 3255 µm). No die, die ID, street etch line or other metal inside. The closest street edge is 55 µm away (lower pair) | Bounding box + 50 µm per side = 184 × 168 µm opening in the ground metal (`tiau_ground_opening.py`, `--margin 50`) |
| How the base complies | Dies on the Mask2 grid (10 mm dies, 10.2 mm pitch, odd 9 × 9 grid with a die at (0, 0); `assemble_wafer_mask2_grid.py`). The wafer-marker rectangles then fall outside every die | Same grid, so every die has the 4 sites at the same local position. `tiau_ground_opening.py` cuts the 4 openings in each source cell |
| Check | `tiau_marker_overlay.py`: every wafer site `outside_dies/clear` (8 / 8) | `tiau_marker_overlay.py`: every chip site `clear` (4 per die). No metal left in any opening, and no die ID touched |

**Bare rectangle in each tone:**
- *Positive*: no polygon inside the rectangle.
- *Etch*: the rectangle is fully exposed.

Check both files after inversion. In particular, the etch scope must expose the wafer-marker rectangles that lie
outside the dies.

**Grid.** A wafer that is not on the Mask2 grid cannot go to NTHU: the fixed markers would land on patterns.
Background and history: [note/tiau_marker_alignment.md](../note/tiau_marker_alignment.md).

### TSRI: A12 marks (pending)

Reference: `Template/TSRI ebeam gds/A12mark_TSRI.gds`, received 2026-11-20. Its geometry, as read from the file:

- **Field**: a 12 × 12 mm field (±6000 µm), layer 0/0.
- **Fine crosses**: four groups of 4 × 4 crosses, 50 µm across with 4 µm lines, on a 100 µm pitch, centred on x, y ∈
  {±3100, …, ±3400} µm.
- **Coarse crosses**: four crosses, 1200 µm long with 4 µm lines, at (±3600, ±4200) µm.
- **Corner marks**: two 500 × 50 µm L marks at (±6000, ±6000) (diagonal), outside a 10 mm die.

To be defined later, before any TSRI delivery:
- metal or opening
- target layer
- per die or per field
- bare-rectangle margin
- whether a wafer-level mark set exists

Then fill this section in with the same detail as NTHU.

## 3. E-beam GDS (per site)

- **One file per site and source:** `<SITE> ebeam gds/<date>_<source>_ebeam_<SITE>.gds`.
- **Junction layers** follow [ebeam_junction.md](ebeam_junction.md): 1–4, or 11–14 when the base sits on 1/0. The
  README and manifest state the layer map, because sites expect different numbers.
- **Base metal**: included as a reference layer only. It must equal the laser-writer positive of the same delivery
  (same source SHA256, vertex-identical after readback).
- **Markers**: the site's marker artwork is included on its own reference layers. For NTHU: wafer markers 3/0 and chip
  markers 2/0, as in Mask2, unless the site asks for other numbers.
- **Junction sizes** come from `code/qubit/jj_resistance_model.py`. The README records the EJ inputs and calibration.

## 4. Folder layout and names

```text
Fabrication output/
  Template/                          tracked: empty site folders + marker references + README
    Laser writer gds/
    NTHU ebeam gds/   TiAu_marker_wafer_4inch.gds               4-inch wafer markers (3/0)
                      TiAu_marker_chip_single_die.gds           chip markers of one die (2/0)
    TSRI ebeam gds/   A12mark_TSRI.gds                          pending
    NCKU ebeam gds/   (no marker)
  <YYYYMMDD>/                        one delivery; local only (git-ignored, files often > 100 MB)
    README.md, manifest.json
    Laser writer gds/                <date>_<source>_positive.gds, _etch.gds (+ _etch.dxf)
    <SITE> ebeam gds/                <date>_<source>_ebeam_<SITE>.gds
```

- **Starting:** copy `Template/` to `<YYYYMMDD>/` using the send date, then delete the unused site folders. Never put
  deliverables in `Template/`.
- **Revisions:** a second delivery on the same date goes to `<YYYYMMDD>_r2`. A delivered folder is never overwritten.
  A replacement goes into a new folder whose README names the one it replaces.
- **`<source>`:** the wafer stem or the Design and revision, e.g. `wafer_summery_mask2grid_open` or `Design3_r11`. No
  spaces.

**Marker references:**

| File | SHA256 |
|---|---|
| `NTHU ebeam gds/TiAu_marker_wafer_4inch.gds` (= Mask2 3/0, 24 polygons) | `7eafe6bec85db93ff614aa0e705673bcb0cdebd43cc3527b3315acddfe989f78` |
| `NTHU ebeam gds/TiAu_marker_chip_single_die.gds` (= Mask2 2/0 within ±5000 µm, 36 polygons) | `d66472dbaec88d81aa89b861ee7e23d2ec0519a2b73eac36d684c76958ab98fb` |
| `TSRI ebeam gds/A12mark_TSRI.gds` | `8d91a98d5b449c86eff7c0abfe6f32dea42ffd85907c6e938ab69a2629149aad` |

The two TiAu files were extracted unchanged from the full photomask `Fluxonium_2Q_bandpass_filter_Mask2.gds`
(SHA256 `9f0ff0793f89429e60a370e76ab3f10285035cee5e390c5a06119eae054e616c`). The combined file is not kept in
`Template/`; use the two split references.

The marker references are fixed: never edit them. A new version is a new file with a new row here.

## 5. README and manifest of a delivery

**`README.md`** (English, short). It covers:
- date, revision, and the folder it replaces
- the e-beam site
- the source design and revision, with its `final_selection` status
- for each file: tone and layer map
- the street grid and die-ID range
- the markers used and their bare-rectangle margin
- checks done, and pending or not-validated items

**`manifest.json`:**

```json
{
  "date": "20261121", "revision": "r1", "replaces": null, "ebeam_site": "NTHU",
  "source": {"gds": "...", "sha256": "...", "design_revision": "...", "final_selection": "pending"},
  "grid": {"chip_um": [10000, 10000], "pitch_um": [10200, 10200], "street_um": 200, "die_ids": "D00-D60"},
  "markers": {"wafer": {"file": "Template/NTHU ebeam gds/TiAu_marker_wafer_4inch.gds", "sha256": "7eafe6be...",
                        "margin_um": 50, "sites_clear": 8},
              "chip": {"file": "Template/NTHU ebeam gds/TiAu_marker_chip_single_die.gds", "sha256": "d66472db...",
                       "margin_um": 50, "sites_clear": 244}},
  "files": [
    {"path": "Laser writer gds/..._positive.gds", "sha256": "...", "tone": "positive", "layers": {"1/0": "metal"}},
    {"path": "Laser writer gds/..._etch.gds", "sha256": "...", "tone": "etch", "scope": "dies-and-streets",
     "made_from": "..._positive.gds"},
    {"path": "NTHU ebeam gds/..._ebeam_NTHU.gds", "sha256": "...", "layers": {"1/0": "base reference",
     "11/0": "small-JJ main", "12/0": "shadow", "13/0": "array main", "14/0": "array shadow",
     "2/0": "TiAu chip markers", "3/0": "TiAu wafer markers"}, "ej_inputs_ghz": {}}
  ],
  "scripts": {}, "checks": {}, "pending": []
}
```

## 6. Checklist before sending

1. **Streets**: present and extended to the wafer edge, and identical in the positive and etch files.
2. **Die IDs**: unique, consecutive from D00, upper-right of each die, correct polarity in both tones.
3. **Tones**: positive and etch come from one source and pass the complement check. The DXF matches the GDS.
4. **Markers**: as the site requires.
   - NTHU: dies are on the Mask2 grid, all 8 wafer sites and every chip site are clear, and each bare rectangle is empty
     in the positive and exposed in the etch.
   - TSRI: not deliverable while pending.
   - NCKU: no markers.
5. **E-beam file**: holds the junction layers and marker references. Its base reference equals the laser-writer
   positive.
6. **Units**: µm with 1 nm precision. The dated folder is new, every file hash is in the manifest, and every
   unconfirmed item is listed under `pending`.
7. **Confirmation**: the source version is user-confirmed per [AGENTS.md](../AGENTS.md). Sending a candidate needs an
   explicit user decision.

Current gap: no script yet assembles a full dated delivery in one step, and the Mask2-grid wafer has no
dies-and-streets etch artwork yet ([note](../note/tiau_marker_alignment.md)). Build each file with the tools listed
above and record it in the manifest.
