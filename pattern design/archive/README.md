# Archived QComponents: not recommended

These parts belong to the old way of building a Purcell filter + readout resonator pair: two separately routed
lines, each with its own TapJunction, joined by a standalone hammer (`HammerCap`) or a finger capacitor
(`InterdigitalCap`). New designs use
[`PurcellReadout`](../resonator/purcell_readout.py) ([method](../../skill/purcell_resonator_design.md),
[layout note](../../note/purcell_readout_parasitics.md)). The files stay importable (`archive/` is on the
`project_paths` search path), so old Design scripts still run.

| File | Class | Replaced by |
|---|---|---|
| `hammer_cap.py` | `HammerCap` | `PurcellReadout` (same heads: `resonator/hammer_geometry.hammer_shapes`) |
| `interdigital_cap.py` | `InterdigitalCap`, `idc_capacitance` | `PurcellReadout` hammer (C_j); `EnclosingCap` (C_p / C_in) |

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
