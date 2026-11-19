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

## E-beam junction overlays on Design3 and Design6 (2026-11-19 r3, in development)

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
