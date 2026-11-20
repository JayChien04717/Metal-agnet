# DesignTestJunction

**Figures:** [fig/](fig/README.md) holds the whole chip ([`chip_r2_fluxonium.png`](fig/chip_r2_fluxonium.png), [`chip_r2_transmon.png`](fig/chip_r2_transmon.png), [`chip_r3_design3.png`](fig/chip_r3_design3.png), [`chip_r3_design6.png`](fig/chip_r3_design6.png)) and one close-up per part, rendered from the current GDS by `code/general/plot_design_figs.py`.

Two standalone process test structures with upper and lower floating probe pads and an
open junction site. The fluxonium style (F) aligns the two facing hook recesses and
reserves room for an array. The transmon style (T) shifts the upper recess across the
stem. The etched letter is at the upper right of each coupon.

| Candidate | Fabrication GDS | Unnotched reference GDS | Layout | Hook detail |
|---|---|---|---|---|
| Fluxonium style | [Fab](gds/DesignTestJunction_2026-11-16_r2_fluxonium_fab.gds) | [Reference](gds/DesignTestJunction_2026-11-16_r2_fluxonium_reference.gds) | Layout (`DesignTestJunction_2026-11-16_r2_fluxonium_layout.png`, removed 2026-11-19; in git history) | Hooks (`DesignTestJunction_2026-11-16_r2_fluxonium_hook.png`, removed 2026-11-19; in git history) |
| Transmon style | [Fab](gds/DesignTestJunction_2026-11-16_r2_transmon_fab.gds) | [Reference](gds/DesignTestJunction_2026-11-16_r2_transmon_reference.gds) | Layout (`DesignTestJunction_2026-11-16_r2_transmon_layout.png`, removed 2026-11-19; in git history) | Hooks (`DesignTestJunction_2026-11-16_r2_transmon_hook.png`, removed 2026-11-19; in git history) |

## Design3 with NCKU-style e-beam junctions (2026-11-20 r4, draft)

This is the current Design3 junction overlay. It takes the Design3 r11 fab base metal unchanged and adds the NCKU EBL
standard components from [`pattern design/ebeam_junction`](../pattern%20design/ebeam_junction/README.md), whose
geometry and naming follow the NCKU final package:

| File | Link |
|---|---|
| Overlay GDS | [GDS](gds/DesignTestJunction_2026-11-20_r4_design3_nckujj_fab.gds) |
| Overlay GDS SHA256 | `d341158238f4b65373c07d96839f34fee5857ac5ce80892a767f9b4303ab559e` |
| Parent | [Design3 r11 fab](../Design3/gds/Design3_2026-11-18_wafer_r11_fab.gds) |
| Sites figure | [figure](notes/fig/DesignTestJunction_2026-11-20_r4_design3_nckujj_sites.png) |
| Record | [JSON](parameters/DesignTestJunction_2026-11-20_r4_design3_nckujj_record.json) |
| Script | [export_ncku_ebeam_junctions.py](layout/export_ncku_ebeam_junctions.py) |

| Site | Component | EJ/h per JJ (input) | Mask widths (finger / cross electrode) |
|---|---|---:|---|
| FX1 (−340, −800) | `FluxoniumJunction`, `side='left'` | 4 GHz | 80 / 92 nm |
| FX2 (340, −800) | `FluxoniumJunction`, `side='right'` | 4 GHz | 80 / 92 nm |
| TC SQUID (0, −1061) | `TransmonJunction`, JJ2 = JJ1 + 30 µm | 11.528 / 11.528 GHz | 80 / 396, 396 nm |

- **Widths**: they come from the NCKU W80 model ([`code/qubit/jj_resistance_model.py`](../code/qubit/jj_resistance_model.py)).
  They are estimates, not measured values.
- **TC EJ (to revalidate)**: the user chose the symmetric NCKU assignment on 2026-11-20. Design3's own CZ model uses the
  asymmetric branches 16.401 / 6.655 GHz, which would give 594 / 199 nm. With the symmetric pair the 9.746 GHz half-flux
  EJ is lost, so the Design3 CZ result needs revalidation.
- **TC contacts**: Design3's ground hooks sit +2 µm from the island hooks, while NCKU Design2's sit +4 µm. The lower
  contact therefore uses `lower_contact_offset='2um'` and `lower_contact_lead_offset='1um'`.
- **Layers**:

  | Layer | Content |
  |---|---|
  | 1/0 | Base metal |
  | 11/0 | Small-JJ main dose (NCKU 1/0) |
  | 12/0 | Shadow (NCKU 2/0) |
  | 13/0 | Array main dose (NCKU 3/0) |
  | 14/0 | Array shadow (NCKU 4/0) |

Checks (geometry only, rerun on the read-back GDS):
- Every hook is filled, and each contact touches only its own electrode.
- Every other shape stays ≥ 1.29 µm (TC) / ≥ 9.19 µm (FX) from the base metal.
- There is no main-dose overlap and no small-JJ main / shadow overlap.
- The base metal is identical after readback.

Not done:
- Design3's corner test junctions get no e-beam JJ.
- No evaporation, HFSS/Q3D or process validation.
- This is not a confirmed release.

```powershell
python -X utf8 DesignTestJunction/layout/export_ncku_ebeam_junctions.py
```

## E-beam junction overlays on Design3 and Design6 (2026-11-19 r3, historical)

The Design3 r3 overlay is superseded by r4 above. These overlays use the archived layout3 components
(`pattern design/archive/*_junction_layout3.py`).

These are draft full-chip GDS files. Each copies the parent fab base metal (1/0) unchanged
and adds `TransmonJunction` (one per SQUID junction) and `FluxoniumJunction` drawings
with the new `hook_contact=True` tapered contacts. Each contact fills its 2 × 4 µm Dolan
hook. The hook mouths are measured on each parent mask. Layers are 11 junction, 12 recess,
13 array, 14 array bar.

| Design | Parent fab GDS | Sites | Overlay GDS | Sites figure | Record |
|---|---|---|---|---|---|
| Design3 | [r11 fab](../Design3/gds/Design3_2026-11-18_wafer_r11_fab.gds) | TC (2 junctions), FX1, FX2 | [GDS](gds/DesignTestJunction_2026-11-19_r3_design3_ebeamjj_fab.gds) | Figure (`DesignTestJunction_2026-11-19_r3_design3_ebeamjj_sites.png`, removed 2026-11-19; in git history) | [JSON](parameters/DesignTestJunction_2026-11-19_r3_design3_ebeamjj_record.json) |
| Design6 | [purcell_resonator_pair_r1 TXT02 fab](../Design6/gds/Design6_2026-11-19_wafer_purcell_resonator_pair_r1_txt02_fab.gds) | Q1, Q2, TC (2 junctions each) | [GDS](gds/DesignTestJunction_2026-11-19_r3_design6_ebeamjj_fab.gds) | Figure (`DesignTestJunction_2026-11-19_r3_design6_ebeamjj_sites.png`, removed 2026-11-19; in git history) | [JSON](parameters/DesignTestJunction_2026-11-19_r3_design6_ebeamjj_record.json) |

Checks (geometry only):
- Every hook is fully covered.
- Each contact touches only its own electrode.
- All other junction shapes stay at least 1.6 µm from the base metal.
- GDS readback confirms the base metal is unchanged.

Not done:
- No EM, fabrication or process validation.
- All junctions use the component default areas. The Design6 asymmetric-SQUID junction sizes are not mapped.
- The Design6 parent is the newest fab mask (`purcell_resonator_pair_candidate`, pending). It is not a user-confirmed release.

```powershell
python -X utf8 DesignTestJunction/layout/export_ebeam_junctions.py --design Design3
python -X utf8 DesignTestJunction/layout/export_ebeam_junctions.py --design Design6
```

The script refuses to overwrite existing outputs. Bump `REVISION` in it for a new export.

See [design.json](design.json) for the candidate index and [layout notes](notes/README.md)
for geometry, validation and fabrication limits. No HFSS, Q3D or quantum simulation was run.
Both r2 mask pairs passed readback of the pad dimensions, hook offsets, etched letters, and fabrication keep-outs. The index records the four GDS hashes. Both fab masks use 5 × 5 µm openings on a 20 µm centre pitch with a 500 µm edge keep-out.

From the repository root, use the `qmetal` environment and a desktop session:

```powershell
python -I -B DesignTestJunction/layout/build.py
```

The [build script](layout/build.py) loads the [geometry inputs](parameters/geometry.json)
and exports both variants. Run [the fluxonium script](layout/test_junction_fluxonium_fab.py)
or [the transmon script](layout/test_junction_transmon_fab.py) with the same Python
command to regenerate one variant. To place a structure in a larger chip, call
`add_test_junction(design, variant, name, fab_contacts=True, pos_x=..., pos_y=..., orientation=...)`
from the build module and export the whole chip. The chip's existing test-junction settings
are independent of these coupons.
