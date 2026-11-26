Only Design-specific QComponents belong here; move mature reusable components to the shared library.

| File | Origin | Use in Design12 |
|---|---|---|
| `interdigital_cap.py` | Copy of `Design8_reflection_purcell_resonator/layout/components/interdigital_cap.py` (2026-11-26, unchanged code) | `IDC_in`, the weak input coupler of the shared λ/4 filter |
| `round_corners.py` | Copy of `Design11_Transmon2Q_Manhattan/layout/components/round_corners.py` (2026-11-26, unchanged code) | R = 5 µm rounding of the Q1 / Q2 / TC core polygons |

`InterdigitalCap` is now used by Design8 and Design12. It is a candidate for promotion to `pattern design/resonator/` (check first for the flat-name clash with the legacy `interdigital_cap` in Design5_6).
