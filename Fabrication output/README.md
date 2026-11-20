# Fabrication output

Every GDS sent for fabrication lives here, one folder per delivery date. The rules are in
[skill/fabrication_output.md](../skill/fabrication_output.md).

Every delivered wafer GDS must have:
- dicing streets
- die IDs
- a positive and an etch (negative) version
- the alignment markers of the e-beam site, each in a bare rectangle

| Folder | Tool / site | Markers |
|---|---|---|
| `Laser writer gds/` | Base metal: `*_positive.gds` and `*_etch.gds` (+ DXF) | Base metal follows the site's marker rule |
| `NTHU ebeam gds/` | NTHU e-beam | TiAu: 4-inch wafer markers + chip markers, bare rectangles (bounding box + 50 µm) |
| `TSRI ebeam gds/` | TSRI e-beam | TSRI A12 marks, **pending** |
| `NCKU ebeam gds/` | NCKU e-beam | None |

- **`Template/`** holds the empty folders and the fixed marker references. Copy it to `<YYYYMMDD>/` to start a delivery,
  and never put deliverables in it.
- **`<YYYYMMDD>/`** (and `<YYYYMMDD>_r2`, …) is one delivery, with its `README.md` and `manifest.json`. These folders
  are local only (git-ignored).
