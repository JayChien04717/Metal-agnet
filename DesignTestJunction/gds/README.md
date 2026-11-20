# GDS layer map

Both file pairs use positive metal on layer 1, µm coordinates and 1 nm precision.
The fabrication mask uses the shared [GDS export convention](../../skill/gds_export.md):
perforated ground on 1/100 and pad metal on 1/10 or 1/11. Each reference contains
unperforated ground and unnotched pads on layer 1. It is a geometry comparison only.
Both file types leave the junction site open.

The 2026-11-19 r3 `*_ebeamjj_fab.gds` files are full-chip overlays: the parent
Design3 / Design6 fab base metal on 1/0, copied unchanged, plus the e-beam junction
drawings on 11/0 (junction and transmon hook contacts), 12/0 (recess), 13/0
(fluxonium array and hook contacts) and 14/0 (array bar).

The 2026-11-20 r4 `DesignTestJunction_2026-11-20_r4_design3_nckujj_fab.gds` uses the NCKU-style components.
Its base metal is the Design3 r11 fab mask on 1/0, unchanged. The NCKU junction layers 1/2/3/4 are written as
11/0 (small-JJ main dose), 12/0 (shadow), 13/0 (array main dose, funnel contacts) and 14/0 (array shadow).
In NCKU's own files the base metal is on 5/0 and the junctions are on 1–4.

See the [artifact index](../README.md) and [geometry checks](../notes/README.md).
