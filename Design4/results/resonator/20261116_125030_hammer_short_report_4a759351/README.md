# Short hammer candidate: arm50 / head30 / gap40 um

Actual full-channel HFSS and local Q3D were rerun. Both main line lengths stay 4.535956670430899 / 4.852581441726434 mm. The requested 50 um arms are measured from each head outer edge to its main-line tap. Head width25 um, CPW16/9 um, ground clearance16 um remain unchanged. Centre spacing is 190 um. Two hammer ground bridges use12 um wire diameter because the17 um landing interval fails the previous15 um wire plus edge-margin check. Other bridges retain their existing definitions. This is simulation geometry, not an adopted fabrication version.

| Quantity | Original | Candidate |
|---|---:|---:|
| Effective J (MHz) | 10.576423 | 12.667352 |
| Bare Purcell frequency (GHz) | 5.980936 | 6.106048 |
| Bare readout frequency (GHz) | 5.951426 | 6.108417 |

Relative J change: +19.77%. These are shared-complex-pole fits of S11/S21 followed by a conditional two-mode inversion, with negligible intrinsic readout loss and real reciprocal coupling. They are not half the visible dip splitting. No physical lengths were changed by fitting. A synthetic known-J recovery check passed. Polynomial background degree sensitivity is 1.31e-07; this does not include EM uncertainty.

Candidate final HFSS Delta S: 0.120970, target0.01. Baseline Delta S: 0.008750. Fine/adaptive history unchanged: True. No direct-frequency or field-mode validation was performed. A5% comparison is only an assistant-selected exploratory screen, not a user-specified acceptance or fabrication criterion.

Local mutual capacitance: 0.488050 fF. Excess shunts: -1.739065 / -1.737005 fF. Signed negative excess values can occur after subtracting a more capacitive straight-CPW reference; they are not standalone physical negative capacitors. Coupling preservation is assessed using full-channel spectra, not by enforcing the old mutual capacitance.

[Geometry and actual spectra](candidate.png) · [Exact results and source runs](comparison.json) · [Run provenance](run.json)
