# DesignTestJunction

Two standalone process test structures with upper and lower floating probe pads and an
open junction site. The fluxonium style (F) aligns the two facing hook recesses and
reserves room for an array. The transmon style (T) shifts the upper recess across the
stem. The etched letter is at the upper right of each coupon.

| Candidate | Fabrication GDS | Unnotched reference GDS | Layout | Hook detail |
|---|---|---|---|---|
| Fluxonium style | [Fab](gds/DesignTestJunction_2026-11-16_r2_fluxonium_fab.gds) | [Reference](gds/DesignTestJunction_2026-11-16_r2_fluxonium_reference.gds) | [Layout](fig/DesignTestJunction_2026-11-16_r2_fluxonium_layout.png) | [Hooks](fig/DesignTestJunction_2026-11-16_r2_fluxonium_hook.png) |
| Transmon style | [Fab](gds/DesignTestJunction_2026-11-16_r2_transmon_fab.gds) | [Reference](gds/DesignTestJunction_2026-11-16_r2_transmon_reference.gds) | [Layout](fig/DesignTestJunction_2026-11-16_r2_transmon_layout.png) | [Hooks](fig/DesignTestJunction_2026-11-16_r2_transmon_hook.png) |

See [design.json](design.json) for the candidate index and [layout notes](notes/README.md)
for geometry, validation and fabrication limits. No HFSS, Q3D or quantum simulation was run.
Both r2 mask pairs passed readback of the pad dimensions, hook offsets, etched letters, and fabrication keep-outs. The index records the four GDS hashes. Both fab masks use 5 × 5 µm openings on a 20 µm centre pitch with a 500 µm edge keep-out.

From the repository root, use the `qmetal` environment and a desktop session:

```powershell
python -I -B DesignTestJunction/layout/build.py
```

The [build script](layout/build.py) loads the [geometry inputs](parameters/geometry.json)
and exports both variants. Run [the fluxonium script](layout/test_junction_fluxonium_fab.py)
or [the transmon script](layout/test_juncion_transmon_fab.py) with the same Python
command to regenerate one variant. To place a structure in a larger chip, call
`add_test_junction(design, variant, name, fab_contacts=True, pos_x=..., pos_y=..., orientation=...)`
from the build module and export the whole chip. The chip's existing test-junction settings
are independent of these coupons.
