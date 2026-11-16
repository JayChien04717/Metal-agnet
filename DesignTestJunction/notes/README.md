# Layout notes

Candidate revision `2026-11-16_r1`. The two variants use the same pads, pocket and landing
stem width. Measurements are in µm in the vertical mask view.

| Geometry | Fluxonium style | Transmon style |
|---|---:|---:|
| Each pad: width (x) × length (y) | 100 × 300 | 100 × 300 |
| Pad corner radius | 10 | 10 |
| Pad body separation | 120 | 120 |
| Landing stem width | 10 | 10 |
| Opening between the uncut facing tips | 30 | 5 |
| Each recess: transverse width × depth | 2 × 4 | 2 × 4 |
| Lower / upper recess centre (x) | 0 / 0 | 0 / +2 |
| Pocket: width (x) × length (y) | 230 × 820 | 230 × 820 |
| Upper-right etched mark | F | T |
| Space reserved for the array | 60 × 10, toward −x | — |

The fluxonium opening follows the existing spacing between the fluxonium landing stems.
The transmon opening follows one branch of the existing TC SQUID, with the upper
electrode taking the role of the original shifted ground-side contact. Both structures
use isolated probe pads. The GDS contains base metal only; the fab supplies the junction
and leads. The openings in the table are mask dimensions; the tunnel barrier geometry is
outside the scope of this layout. The fluxonium array envelope extends 60 µm to
the left of the junction centre with 55 µm clearance to ground, and does not overlap
either stem. The array and SQUID loop are omitted from the base metal. See the
[parameter sources](../parameters/provenance.json) for the reference functions and
design choices.

The standalone 1.5 × 2.0 mm frame gives at least 590 µm from pocket to die edge and the
pocket gives 50 µm minimum ground-to-pad clearance. The files use µm coordinates and
1 nm precision.
Fabrication ground carries the lab's hole pattern and pocket keep-out. The build fills
holes clipped by the pocket keep-out and any hole that crosses the 500 µm die-edge
keep-out. See the
[layer map](../gds/README.md). The unnotched references are layout comparisons only.

The build reads the GDS files back and checks their outline, layer, resolution, pad
dimensions, connected ground, two floating pads, open junction centre and the exact two
recess cuts. It checks both etched letters and the fluxonium array envelope. All retained
holes are 5 × 5 µm, at least 50 µm from the pocket and letter, and at least 500 µm from
the die edge. The layout used the lab
`qmetal` environment (Qiskit Metal 0.7.4, gdstk 1.0.0, gdspy 1.6.13,
Shapely 2.1.2). Laboratory process rules for the separate junction layers and the eventual
chip placement still need validation. No HFSS or Q3D run applies to this process layout.
