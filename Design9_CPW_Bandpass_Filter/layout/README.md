# Design9_CPW_Bandpass_Filter Layout

- Entry script: [`build.py`](build.py)
- Uses shared QComponents from `pattern design/`: `BondPad` (`bondpad.py`), `TapJunction` (`resonator/tap_junction.py`), and standard `qiskit_metal` `ShortToGround`, `RouteStraight`, `RouteAnchors`.
- All dimensions are loaded from [`../parameters/geometry.json`](../parameters/geometry.json).
