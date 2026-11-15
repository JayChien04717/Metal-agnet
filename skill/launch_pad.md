# Design6 wire-bond launch pad standard

Use the shared QComponent in [bondpad.py](../pattern%20design/bondpad.py). Use
`add_bond_pad` for a nominal PCB anchor. It offsets the taper origin inward along the
pad normal. Use `BondPad` directly for fixtures where the taper origin is specified explicitly.

| Parameter | Design6 standard |
|---|---:|
| Bonding metal width | 400 µm |
| Bonding metal length along the route | 300 µm |
| Side clearance from metal to ground | 200 µm on each side |
| Clearance behind the pad | 50 µm |
| Linear taper length | 300 µm |
| Straight lead between taper and `tie` pin | 30 µm |
| Inward shift from nominal PCB anchor | 250 µm |

The `tie` pin uses the connected CPW width and gap. The `bond` pin is on the
back edge of the 400 µm pad and faces the chip edge across the 50 µm gap.
Design6 uses 10/6 µm for control lines and 16/9 µm for readout. Other layouts
carry their existing CPW widths and gaps through to the launch.
Preserve each design’s assigned ports, orientation and net connectivity.
On a 10 mm chip with nominal anchors at ±4.5 mm, the etch stops 100 µm from the
die edge. Design2 retains its GDS reference port coordinates along the edge.
Its 9.9 mm metal field stops 50 µm inside the 10 mm die. The Design4 single-channel
fixture uses explicitly positioned taper origins and keeps its 5.2 × 4.3 mm outline.

The pad taper uses the Design6 300 µm geometry and is exempt from the generic taper-angle check.
Inspect the pad and both CPW slots in the exported GDS. Preserve existing solver results with their
original launch geometry. Re-render the model for any new whole-chip or driven response validation.
