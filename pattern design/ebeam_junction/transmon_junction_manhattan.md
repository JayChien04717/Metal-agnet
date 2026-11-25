# `TransmonJunctionManhattan`

Two-layer e-beam mask drawing for a Manhattan-style Josephson junction (single JJ or 2-JJ SQUID) on a transmon or tunable coupler. Geometry matches the 45°-rotated-pad Manhattan standard (4 × 8 µm pads, orthogonal x/y fingers crossing at `(+2, -2) µm` with 1 µm stubs, 5 × 9 µm pad undercut/bandage boxes, and 1 µm finger-tip undercut boxes).

## Orientation rule

All Manhattan junctions on a chip must keep their fingers strictly along the chip `x` and `y` axes so one two-angle shadow evaporation covers every junction. When `squid='True'`, branch 2 (`jj2`) is branch 1 (`jj1`) translated along `+x` by `squid_spacing` — never mirrored.

## Options

| Option | Default | Role |
|---|---|---|
| `jj1_finger_x_width` | `'190nm'` | Horizontal finger width (y-extent) of JJ1 |
| `jj1_finger_y_width` | `'170nm'` | Vertical finger width (x-extent) of JJ1 |
| `jj2_finger_x_width` | `'190nm'` | Horizontal finger width of JJ2 (when `squid='True'`) |
| `jj2_finger_y_width` | `'170nm'` | Vertical finger width of JJ2 (when `squid='True'`) |
| `finger_stub_length` | `'1.0um'` | Finger extension past the crossing region |
| `cross_x`, `cross_y` | `'2.0um'`, `'-2.0um'` | Fixed centre of the finger crossing in the JJ local frame |
| `pad_length`, `pad_width` | `'8.0um'`, `'4.0um'` | Size of each 45°-rotated landing pad |
| `pad_inner_gap` | `'4.456um'` | Inner separation between the two 45°-rotated pads along the +45° axis |
| `pad_undercut_margin` | `'0.5um'` | Margin around each pad for the undercut/bandage rectangle (→ 5 × 9 µm) |
| `tip_undercut_length` | `'1.0um'` | Length of each finger-tip undercut box |
| `tip_undercut_overlap` | `'0.1um'` | Overlap of the tip undercut box onto the finger tip |
| `tip_undercut_margin` | `'0.125um'` | Side margin added to each side of the finger for the tip undercut box |
| `squid` | `'True'` | `'True'` draws two translated JJs (`jj1`, `jj2`); `'False'` draws one centred JJ |
| `squid_spacing` | `'30um'` | Centre-to-centre x separation between `jj1` and `jj2` |
| `layer_jj` | `'1'` | GDS layer for the main Manhattan JJ polygon (use e.g. `'11'` or `'13'` on a chip) |
| `layer_jj_undercut` | `'2'` | GDS layer for pad undercut/bandage and tip undercut boxes (use e.g. `'12'` or `'18'` on a chip) |
| `jj_metal` | `'True'` | `'False'` omits all polygons |

## Pins

For each drawn branch `N` (`1`, and `2` when `squid='True'`):
- `top_jjN`: centre of the upper-right landing pad, outward normal along +45°.
- `bottom_jjN`: centre of the lower-left landing pad, outward normal along −135°.
