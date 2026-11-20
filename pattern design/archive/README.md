# Archived QComponents: not recommended

The two `*_junction_layout3.py` e-beam drawings are explained below the table. The capacitor parts belong to the old way of building a Purcell filter + readout resonator pair: two separately routed
lines, each with its own TapJunction, joined by a standalone hammer (`HammerCap`) or a finger capacitor
(`InterdigitalCap`). New designs use
[`PurcellReadout`](../resonator/purcell_readout.py) ([method](../../skill/purcell_resonator_design.md),
[method note section 6](../../note/purcell_filter_method.md#6-purcellreadout-layout-hammer-only-coupling)). The files stay importable (`archive/` is on the
`project_paths` search path), so old Design scripts still run.

| File | Class | Replaced by |
|---|---|---|
| `hammer_cap.py` | `HammerCap` | `PurcellReadout` (same heads: `resonator/hammer_geometry.hammer_shapes`) |
| `interdigital_cap.py` | `InterdigitalCap`, `idc_capacitance` | `PurcellReadout` hammer (C_j); `EnclosingCap` (C_p / C_in) |
| `fluxonium_junction_layout3.py` ([doc](fluxonium_junction_layout3.md)) | `FluxoniumJunction` (layout3 drawing) | NCKU-style `ebeam_junction/fluxonium_junction.py` ([README](../ebeam_junction/README.md)) |
| `transmon_junction_layout3.py` ([doc](transmon_junction_layout3.md)) | `TransmonJunction` (mirrored SQUID) | NCKU-style `ebeam_junction/transmon_junction.py` |

The two `*_junction_layout3.py` modules were moved here on 2026-11-20 with unchanged bytes, so the SHA256 values recorded
by Design3 r11 and DesignTestJunction r3 still match. They use the old option names (`junction_gap`, `junction_width`,
`finger_width`, `mirror`, `hook_contact`) and a mirrored SQUID, which the NCKU standard replaces with a translated pair
in which every finger is on the left. These scripts now import them from `archive/`:

- `Design3/layout/add_ebeam_junctions.py`
- `DesignTestJunction/layout/export_ebeam_junctions.py`
- `code/general/check_fluxonium_junction.py`
- `code/general/check_transmon_junction.py`
- `code/general/plot_junction_parameters.py`

Only their import lines changed. The new mask is built with the NCKU-style components.

## Why not

In the old layouts the filter and resonator ran side by side for 1–1.6 mm at 190–250 µm, and the stub + pocket
etch cut the ground between them. The Q3D hammer C_j then explained only about a quarter of the HFSS J:
- Design4 at 200 µm: hammer closed form 2.9 MHz against HFSS 12.0 MHz.
- Design6: 0.66 fF in Q3D against about 2.2 fF equivalent.

With `PurcellReadout`, the Design8 driven S21 gives J within 1 % of the Q3D hammer value
([Design8](../../Design8/README.md)).

## If you still use them

The capacitance matrix of the hammer coupon is not enough. Correct the model with all of the following:
1. Coupled-line sections for every side-by-side run. Take c′_m from a two-line Q3D or
   `code/resonator/cpw_pair_cap2d.py`, and include the inductive part (k_L ≈ k_C).
2. The rest of the full-channel Q3D filter–resonator mutual, including tees and bends. The local coupon is not
   enough.
3. Ground straps across both stubs and at every T, in Q3D and in HFSS. Without them, the S21 fit puts the
   ground-return coupling into C_j (J was 3× too large in Design4 at 250 µm).
4. HFSS for J and f_p. Treat the circuit model as a trend only (skill/purcell_resonator_design.md, legacy section).
