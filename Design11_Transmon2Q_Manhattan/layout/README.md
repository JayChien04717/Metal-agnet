# Layout folder (`Design<N>/layout/`)

Follow [skill/layout_structure.md](../../../skill/layout_structure.md):

* `placement.json` — declarative component placement (qubits, couplers, tees, launch pads, site markers such as `TSRIA12Mark` / `TiAuChipMarker`, coupons, die ID).
* `routes.json` — CPW routes, explicit anchors, target filleted resonator lengths, and meander rules.
* `build.py` — single entry point `build_design(variant='sim'|'process'|'fab', overrides=None)` reading `placement.json` and `routes.json`, solving resonator lengths via `code/layout/route_length_solver.py`, and verifying keep-outs via `code/layout/check_keepouts.py`. Never import simulation modules here.
* `components/` — optional Design-specific `QComponent` subclasses only; import shared components from root `pattern design/`.
