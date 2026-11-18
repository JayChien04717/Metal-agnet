# Standalone junction array GDS

[fluxonium_junction.gds](fluxonium_junction.gds) is the requested standalone
multilayer junction/array drawing, exported from
[FluxoniumJunction](../pattern%20design/fluxonium/fluxonium_junction.py).
It is a draft for inspection, without a chip ground plane, cheesing or corner
test coupons. The full parameter snapshot, source hash and read-back geometry
checks are in [export_record.json](export_record.json). No fabrication release
or electrical performance is validated; simulations have not been rerun.

[export_gds.py](export_gds.py) reproduces the export in the `qmetal`
environment using the shared `export_component_polygons` helper in
[gds_tools.py](../code/general/gds_tools.py). It refuses to overwrite an existing
output; choose a new filename before re-exporting a changed candidate.

Two additional draft comparison exports are available:

- [Array count variation](fluxonium_junction_n30_gap075nm.gds), with its
  [parameter snapshot and checks](fluxonium_junction_n30_gap075nm_record.json).
- [Shared junction gap variation](fluxonium_junction_n43_gap100nm.gds), with its
  [parameter snapshot and checks](fluxonium_junction_n43_gap100nm_record.json).

These are proposed comparison inputs, not final confirmed versions. The CLI
accepts `--num-cells` (applies to both arms), `--junction-gap` (with units),
`--junction-width` (right-hand junction bar y width, with units) and
`--output` (new file path). The three existing exports predate
`junction_width` and use its 0.25 um default geometry; their recorded
`component_sha256` refers to the earlier component source. The shared arrow clearance parameter in the
QComponent is `junction_gap`; the per-arm rectangle count is `num_cells`.
