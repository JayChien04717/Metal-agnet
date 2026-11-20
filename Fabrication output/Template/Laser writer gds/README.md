# Laser writer

This folder holds base metal only, from one positive source:
- `<date>_<source>_positive.gds`: metal = polygons
- `<date>_<source>_etch.gds`: exposed = polygons, made with `code/general/invert_wafer_tone.py`
- `<date>_<source>_etch.dxf`: when the writer needs DXF

Both tones carry the dicing streets, the die IDs and the e-beam site's marker bare rectangles: empty in the positive,
exposed in the etch. There are no junction layers. Rules:
[skill/fabrication_output.md](../../../skill/fabrication_output.md#1-required-content-of-every-delivered-wafer-gds).
