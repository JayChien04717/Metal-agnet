# E-beam junction masks: EJ → resistance → area → drawn pattern

This skill covers drawing the Dolan-bridge junction layers for a fluxonium (small JJ + JJ array) or a transmon / coupler
SQUID, sizing them from EJ, and placing them on a chip mask. Everything needed is in three files:

| File | Role |
|---|---|
| [`code/qubit/jj_resistance_model.py`](../code/qubit/jj_resistance_model.py) | EJ ↔ Rn ↔ overlap area ↔ mask widths (standard library only) |
| [`pattern design/ebeam_junction/fluxonium_junction_dolan.py`](../pattern%20design/ebeam_junction/fluxonium_junction_dolan.py) | `FluxoniumJunctionDolan`: small JJ + array + funnel contacts ([options](../pattern%20design/ebeam_junction/fluxonium_junction_dolan.md)) |
| [`pattern design/ebeam_junction/transmon_junction_dolan.py`](../pattern%20design/ebeam_junction/transmon_junction_dolan.py) | `TransmonJunctionDolan`: SQUID or single JJ ([options](../pattern%20design/ebeam_junction/transmon_junction_dolan.md)) |

The checks are in [`code/ebeam/check_ebeam_junctions.py`](../code/ebeam/check_ebeam_junctions.py) and the parameter
figures come from [`code/ebeam/plot_ebeam_junction_parameters.py`](../code/ebeam/plot_ebeam_junction_parameters.py).
The geometry and parameter names were adopted from the NCKU fabrication team's EBL standard (2026-11). The old
layout3 drawings are in `pattern design/archive/` and are not used for new masks.

## 1. Geometry rules

1. **Orientation.** Every small JJ has its vertical finger on the left and its horizontal cross electrode on the right,
   so all junctions share one evaporation direction.
   - In a SQUID, branch 2 is branch 1 translated by `squid_spacing`; it is never mirrored.
   - The fluxonium `side='right'` drawing is rotated by 180° and re-routed so that its finger stays on the left. It is
     not a mirror image.
2. **Fixed edges.** When widths change, the finger's right edge and the cross electrode's top edge stay fixed.
   - The electrode's lower edge sets `jj_cross_electrode_width`.
   - The finger grows to the left.
3. **Shadow (undercut) openings.**
   - They keep `shadow_clearance` = 75 nm from the cross electrode.
   - Small-JJ main dose never overlaps shadow.
   - Small-JJ and array main doses only touch at shared edges.
4. **Bridge gaps.**
   - `jj_bridge_gap` is 248 nm (fluxonium) and 249 nm (transmon).
   - The array gap `array_bridge_gap` is 0.202 µm.
5. **Contacts.**
   - Each contact body covers its base-metal hook or bay with a 1 µm margin on each side.
   - The body runs 1 µm past the hook bottom.
   - Use 4 µm for a 2 µm rectangular hook and 7.6 µm for a 5.6 µm rounded bay.
   - The fluxonium funnel mouth (`array_contact_lead_width`) is 2 µm.
6. **Layers** (component defaults):

   | Layer | Content |
   |---|---|
   | 1 | Small-JJ main dose |
   | 2 | Shadow |
   | 3 | Array main dose (cells, return, funnel contacts) |
   | 4 | Array shadow |

   - The components never draw base metal.
   - On a chip whose base is on 1/0, write them on 11–14 instead (`layer_jj`, `layer_jj_shadow`, `layer_array`,
     `layer_array_shadow`).
7. **Scope.** All values are planar mask dimensions. The junction area after evaporation is a model estimate (Section 2).

## 2. EJ → resistance → area → mask width

Every quantity is per single JJ, including each branch of a SQUID. EJ is a design input (user or quantum model), never an
HFSS result.

```text
Rn [kOhm] = 140.433954 / (EJ/h [GHz])                Ambegaokar-Baratoff, Al, Delta = 180 ueV, T -> 0
A_actual  = RA / Rn,   RA = 625 Ohm um^2             one measured junction: 0.100 x 0.250 um actual, about 25 kOhm
A_actual  = (jj_finger_width + 30 nm) * (jj_cross_electrode_width + 70 nm)
=>  A_actual [um^2] = 0.0044505 * EJ/h [GHz]
```

**Which parameter sets the area.**
- `jj_finger_width`, the vertical finger, is fixed at **80 nm** for every small JJ. It is predicted to fabricate at
  110 nm.
- The area is set only by **`jj_cross_electrode_width`**, the y width of the horizontal electrode. On `TransmonJunctionDolan`
  it is set per branch: `jj1_cross_electrode_width`, `jj2_cross_electrode_width`.

```text
jj_cross_electrode_width = A_actual / 110 nm - 70 nm      (rounded to the 1 nm grid)
```

"Actual" widths are predicted fabricated sizes; never draw them. For an EJ that would need more than the drawing limit
(Section 5), widen the finger with `finger_from_ej(ej, cross_width)`.

### Area vs resistance (finger 80 nm)

| `jj_cross_electrode_width` (nm) | Actual area (µm²) | Rn (kΩ) | EJ/h (GHz) |
|---:|---:|---:|---:|
| 80 | 0.01650 | 37.879 | 3.707 |
| 100 | 0.01870 | 33.422 | 4.202 |
| 120 | 0.02090 | 29.904 | 4.696 |
| 150 | 0.02420 | 25.826 | 5.438 |
| 200 | 0.02970 | 21.044 | 6.673 |
| 250 | 0.03520 | 17.756 | 7.909 |
| 300 | 0.04070 | 15.356 | 9.145 |
| 350 | 0.04620 | 13.528 | 10.381 |
| 400 | 0.05170 | 12.089 | 11.617 |
| 450 | 0.05720 | 10.927 | 12.853 |
| 500 | 0.06270 | 9.968 | 14.088 |
| 600 | 0.07370 | 8.480 | 16.560 |

To print the same table for another finger width, pass `--finger` (example below uses 100 nm):

```bash
python -I code/qubit/jj_resistance_model.py --sweep --finger 100
```

### SQUID from its asymmetry

For an ideal SQUID (loop inductance neglected):

```text
EJ(Phi)  = EJ_sum * sqrt(cos^2(pi Phi/Phi0) + d^2 sin^2(pi Phi/Phi0))
EJ_sum   = EJ1 + EJ2 = EJ(0),     d = (EJ1 - EJ2) / EJ_sum = EJ(Phi0/2) / EJ(0)
EJ1, EJ2 = EJ_sum (1 + d) / 2,  EJ_sum (1 - d) / 2,       A1 / A2 = (1 + d) / (1 - d)
```

Example: EJ(0) = 23.056 GHz and EJ(Φ0/2) = 9.746 GHz give d = 0.4227, A1 / A2 = 2.465, and these branches:

| Branch | EJ/h (GHz) | Rn (kΩ) | Area (µm²) | `jjN_cross_electrode_width` |
|---|---:|---:|---:|---:|
| JJ1 | 16.401 | 8.563 | 0.07299 | 594 nm |
| JJ2 | 6.655 | 21.102 | 0.02962 | 199 nm |

Size each branch from its own EJ, never from an idle or CZ effective EJ. A symmetric SQUID (d = 0) has zero EJ at half
flux. Which side gets the larger junction is a layout choice.

### Calibration status

RA comes from one measured junction, and the +30 nm / +70 nm growth values are process assumptions, so every size here
is an estimate. Room-temperature probe resistances need a process calibration before they can be compared with the
low-temperature Rn. After fabrication, fit the model to the test-key resistances and pass the fitted values to every
function: `ra`, `finger_growth_nm`, `cross_growth_nm`, `delta_uev`. Record the fit in the Design.

## 3. Workflow

1. Take the per-JJ EJ from the Design (provenance per [parameter_provenance.md](parameter_provenance.md)).
2. Compute the options:

   ```python
   import sys; sys.path[:0] = ['code/qubit', 'pattern design/ebeam_junction']
   import jj_resistance_model as jjm
   from fluxonium_junction_dolan import FluxoniumJunctionDolan
   from transmon_junction_dolan import TransmonJunctionDolan

   fx = jjm.size_from_ej(4.0).options()                       # {'jj_finger_width': '80nm', 'jj_cross_electrode_width': '92nm'}
   tc = jjm.squid_sizes(23.056, 0.4227)['options']            # or jjm.squid_options(ej1, ej2)
   ```

   The same numbers from the command line:

   ```bash
   python -I code/qubit/jj_resistance_model.py --ej 4 --squid-flux 23.056 9.746
   ```

3. **Measure the hooks** on the parent fab mask's base metal: hook x, the mouth y, and the hook depth and width. Then
   place the components:
   - **Fluxonium**:
     - `pos` = the contact centre (the array rows sit at y = ±4 µm).
     - `side='left'` puts the array towards −x and `side='right'` towards +x.
     - `array_contact_end` = hook mouth + hook depth + 1 µm (20 µm for a 15 µm mouth with 4 µm hooks).
   - **SQUID**:
     - `pos` = the SQUID centre, with y = 0 midway between the island-hook mouth (+2.5 µm) and the ground-hook mouth
       (−2.5 µm).
     - The upper patch is centred on the island hook. Set `lower_contact_offset` to the ground hook's x offset.
     - If that offset is below about 3 µm, also set `lower_contact_lead_offset` so the stem stays ≥ 75 nm from the
       shadow openings. For example, for a +2 µm hook use `'2um'` / `'1um'`. The component refuses a stem that is too
       close.
   - **Single test junction**: `squid=False`. Use `upper_contact_offset` when the upper bay is offset (a 0.5 µm ledge
     joins the finger).
4. **Export:** follow the fab GDS workflow in [gds_export.md §2.1](gds_export.md#21-fab-gds-workflow).
   - Place junctions inside `build_design()` using [`code/ebeam/ebeam_junction_placement.py`](../code/ebeam/ebeam_junction_placement.py) `add_ebeam_junctions()`, which reads `contact_hooks()` from the fab components, and export with `export_gds(..., fab=True, datatype0_layers=(11, 12, 13, 14))`.
   - Verify every hook is filled by its own contact only, all other shapes stay off base metal, and base metal is unchanged.
   - Legacy standalone GDS overlay scripts are historical only and not used for new exports.
5. **Run the checks**:

   ```bash
   python -I -X utf8 code/ebeam/check_ebeam_junctions.py
   ```

   Re-render the parameter figures if the geometry changed:

   ```bash
   python -I -X utf8 code/ebeam/plot_ebeam_junction_parameters.py
   ```

## 4. What the components are not

- They are mask drawings only. They contain no lumped junction, inductance or EJ / EL model and do not replace
  `FluxoniumArrayJJDolan` or the HFSS lumped-L models.
- `num_array` counts drawn cells per arm, not an independently established number of electrical junctions.
- Chip masks stay junction-free by default (`jj_metal=False`, [gds_export.md](gds_export.md)). Adding junction layers to
  a fabrication export is an explicit Design choice.
- No evaporation simulation, HFSS / Q3D run with the drawn junctions, or process qualification is implied.

## 5. Limits

These are geometric limits of the drawing, not process rules. The components raise `ValueError` beyond them.

| Limit | Value |
|---|---|
| `jj_cross_electrode_width`, finger must span the electrode | < 0.687 µm (about EJ/h < 18.7 GHz at 80 nm) |
| `jj_cross_electrode_width`, shadow floor | < 0.862 µm |
| `array_row_width` | ≤ 2 µm |
| Fluxonium `array_contact_body_width` | Must cover the funnel mouth |
| SQUID `squid_spacing` | The two branches must not overlap |
| Upper contact body | Must reach the finger's right edge |
