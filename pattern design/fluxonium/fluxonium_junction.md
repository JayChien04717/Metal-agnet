# Fluxonium junction mask component

`FluxoniumJunction` in [fluxonium_junction.py](fluxonium_junction.py) is a
standalone QComponent for the multilayer junction/array drawing requested on
2026-11-18. It does not modify any existing chip assembly or adopted result.
This is a draft geometry candidate; full parameter and fabrication confirmation
are pending. No HFSS, Q3D, quantum-model or fabrication validation has been run.

## Source and intentional changes

The user supplied `C:/Users/QEL/Downloads/layout3 ftf.gds` and explicitly asked
to reference it. SHA-256:
`6546a6202ce78dc716b100501d88884f8ed6f5854e17b0e59d962bce10ea2b26`.
The source unit is 1 um, precision 1 nm. The selected device is the left
fluxonium at (-485, -1000) um in cell `TOP`. Complete polygons on layers 1-4
contained within x=(-545,-475), y=(-1017,-983) um are retained, translated by
(485,1000) um, in [the reference snapshot](fluxonium_junction_reference.json).
The snapshot was missing from the working tree on 2026-11-18 and was rebuilt
from the same source file (SHA-256 verified identical) with that procedure.
The snapshot makes component construction independent of the Downloads file.
Source fabrication process, electrical layer interpretation and convergence
evidence are unknown. Geometry is evidence of mask coordinates only.

The source has 43 rectangular cells per arm, at approximately 1.2 um pitch,
with 2 um transverse width. Their longitudinal widths range from 0.998 to
1.000 um. The new component deliberately uses one fixed 0.998 um width for
every cell, matching the user's marked example. This uniformization is not an
exact reproduction of the source's nanometre variations. The source's main
junction has unequal clearances (including about 0.050, 0.074, 0.075 and
0.076 um). All relevant clearances around the finger and the right bridge are
instead controlled by one `junction_gap`; its 0.075 um default is a proposed
input based on the source, not a separately confirmed fabrication rule.
The red main-junction polygons are retained exactly at the default
`junction_width`; the adjacent purple polygons are adjusted.

`junction_width` (default 0.25 um, the source value) sets the y width of the
right-hand horizontal red junction bar, resized symmetrically about its source
centre line y = 2.571 um; its x extent (3.732-5.73 um) is fixed. The bridge
recesses follow its edges with `junction_gap`, and the two recess strips beside
the bar near the right contact keep their fixed source 0.075 um clearance.
This width enters the junction overlap area, but the actual area depends on the
evaporation process, which is not modelled here. Its default is a proposed
input from the source, not a confirmed fabrication target.

All rebuilt boxes are checked for strictly positive size, so an infeasible
`junction_gap`/`junction_width` raises `ValueError` rather than silently
producing inverted shapes (the former `junction_gap < 0.25 um` check allowed
gaps above ~0.124 um that inverted the recess between finger and bar).
Geometric feasibility is not a process rule. Original layers 1, 2, 3, 4 are preserved by default.

## Usage

Add this directory to `sys.path` in the consuming layout script, then:

```python
from qiskit_metal import designs
from fluxonium_junction import FluxoniumJunction

design = designs.DesignPlanar()
device = FluxoniumJunction(design, 'FJ', options=dict(
    num_cells=43,             # PER ARM: 86 rectangles total
    cell_length='0.998um',    # along the array, x
    cell_width='2um',         # transverse, y
    cell_gap='0.202um',       # pitch = cell_length + cell_gap
    junction_gap='0.075um',   # shared clearances indicated by the user
    junction_width='0.25um',  # y width of the right-hand junction bar
    pos_x='0um', pos_y='0um', orientation=0,
))
```

For quantity changes, edit only `num_cells`. Cells are never stretched to fit
a fixed total length. The right contacts and main junction stay fixed; the
left return moves by the change in array length. Array bar lengths change to
cover the cells; the bar transverse excess remains 0.05 um on each side.
The return retains a 0.202 um clearance from the last cell. The origin is the
source device's contact center, not the changing overall bounding-box center.
Pins `top` and `bottom` are at local (0, +/-16) um and face outwards.

Optional layer mapping: `layer`, `layer_recess`, `layer_array`,
`layer_array_bar`. These are original mask drawing layers; their overlap must
not be interpreted as galvanic connection without process information. The
component adds positive polygons without a ground plane or etched pocket.
It contains no lumped-junction sheets, inductance, EJ or EL model. Do not use
it as a replacement for `FluxoniumArrayJJ` in an HFSS calculation.
`num_cells` counts drawn cells, not an independently established number of
electrical junctions. `jj_metal=False` omits the entire junction component,
consistent with the repository's default junction-free chip masks. Keeping
this junction drawing in a fabrication export requires an explicit Design
choice, plus validation of the process and layer mapping.

Aliases `fluxonium_junction` and `fluxonium_juncion` are also available.

## Validation and reproduction

Use the documented `qmetal` Python:

```powershell
& 'C:/Users/QEL/miniforge-pypy3/envs/qmetal/python.exe' code/general/check_fluxonium_junction.py
```

The check covers 12 combinations of count (1,12,43,80) and shared gap
(0.05,0.075,0.1 um), 9 combinations of `junction_width` (0.15,0.25,0.35 um)
and gap, positive valid polygons, invariant cell dimensions, six
measured common clearances, unchanged red source polygons, transformed pins,
disabled geometry, and rejection of invalid dimensions/counts. It passed on
2026-11-18 (rerun after adding `junction_width`). This does not establish
fabrication connectivity, minimum process
features, evaporation directions, junction areas or electrical performance.
Nominal pin width is the source's 1.002 um contact width; actual port process
and integration into any chip remain unvalidated.

The same check accepts `--preview-dir <directory>` to save a draft
`MetalGUI.screenshot` of the junction geometry (`jj_metal=True`). This preview
is specifically a junction-mask drawing, not a junction-free chip mask or a
fabrication deliverable. The separately requested standalone draft GDS and its
read-back checks are indexed in [test_array](../../test_array/README.md).
