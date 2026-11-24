# TSRI e-beam: A12 marks

`A12mark_TSRI.gds` is TSRI's mark artwork: a 12 × 12 mm field on layer 0/0. It is reproduced by `a12_groups()` in
[pattern design/marker/tsri_a12_mark.py](../../../pattern%20design/marker/tsri_a12_mark.py) (checked by
`code/general/check_marker_references.py`). Never edit this reference.

Placement plan (2026-11-24, not yet confirmed by TSRI): the fine-cross groups and the four coarse crosses inside every
die (`TSRIA12Mark`, layer `2/0` of the e-beam GDS) with a bare rectangle (bounding box + 50 µm) in the base metal at
every group; the two corner L marks belong to the dicing streets and are drawn with the street artwork. Until TSRI
confirms, TSRI deliveries are review drafts. Rules and procedure:
[skill/fabrication_output.md](../../../skill/fabrication_output.md#tsri-a12-marks).
