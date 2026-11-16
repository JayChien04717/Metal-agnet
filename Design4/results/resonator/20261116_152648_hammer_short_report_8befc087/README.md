# Short hammer candidate: arm50 / head30 / width25 / gap40 um

Actual full-channel HFSS and local Q3D were rerun. Both main line lengths stay 4.535956670430899 / 4.852581441726434 mm. Each 50 um arm is measured from its head outer edge to the main-line tap. Head width is 25 um; CPW16/9 um and ground clearance16 um remain unchanged. Centre spacing is 190 um. The arm follows (centre spacing - gap - 2*head width)/2. Two local ground bridges use12 um wire diameter with landing and etch checks recorded in the candidate geometry audit. This is simulation geometry, not an adopted fabrication version.

| Quantity | Original | Candidate |
|---|---:|---:|
| Effective J (MHz) | 10.576423 | 12.752274 |
| Bare Purcell frequency (GHz) | 5.980936 | 6.107167 |
| Bare readout frequency (GHz) | 5.951426 | 6.110568 |

Relative J change: +20.57%. These are shared-complex-pole fits of S11/S21 followed by a conditional two-mode inversion, with negligible intrinsic readout loss and real reciprocal coupling. They are not half the visible dip splitting. No physical lengths were changed by fitting. A synthetic known-J recovery check passed. Polynomial background degree sensitivity is 1.19e-07; this does not include EM uncertainty.

Candidate final HFSS Delta S: 0.005202, target0.01. Baseline Delta S: 0.008750. Fine/adaptive history unchanged: True. No direct-frequency or field-mode validation was performed. A5% comparison is only an assistant-selected exploratory screen, not a user-specified acceptance or fabrication criterion.

Local mutual capacitance: 0.488050 fF. Excess shunts: -1.739065 / -1.737005 fF. Signed negative excess values can occur after subtracting a more capacitive straight-CPW reference; they are not standalone physical negative capacitors. Coupling preservation is assessed using full-channel spectra, not by enforcing the old mutual capacitance.

[Geometry and actual spectra](candidate.png) · [Exact results and source runs](comparison.json) · [Run provenance](run.json)
