# Actual 400 / 800 um hammer-arm HFSS S21

Full Design4 hammer channels were newly simulated with each head-outer-edge to main-line tap arm set to 400 or 800 um. Main filter/resonator lengths remained 4.535956670430899 / 4.852581441726434 mm. Hammer gap/height/width remained 42 / 65 / 25 um; CPW remained 16 / 9 um.

The straight-arm change increases main-line centre spacing to 0.892 / 1.692 mm. This also changes distributed inter-line coupling. The same grounding policy is retained (30 automatic bonds and 9 explicit straps); this is simulation geometry, not a fabrication release.

Raw full-arm model: the entire corresponding Q3D coupon is a pi network at the two tap planes, without additional arm transmission lines. Deembedded model: subtract matched CPW reference and represent the actual tee-to-hammer-centre arms as transmission lines (0.446 / 0.846 mm each). Both use the same fixed main lengths and unfitted Cin/Cp extraction. Model S21 includes the deterministic renderer port-polarity conversion.

| Actual arm per side (um) | HFSS dips (GHz) | Raw-coupon dips (GHz) | Deembedded-model dips (GHz) | Last Delta S |
|---|---|---|---|---|
| 400 | 5.662875, 5.727450 | 5.721450, 5.736625 | 5.724725, 5.740100 | 0.019345 |
| 800 | 5.323525, 5.420950 | 5.404825, 5.431725 | 5.394425, 5.422025 | 0.022242 |

Requested HFSS Delta S is 0.01. Pass count is not convergence. Fine sweeps reuse each case's adaptive setup; different geometries do not share a mesh. Q3D/HFSS cross-solver mesh reuse is not claimed. The 400/800 um Q3D extractions and unchanged Cin/Cp coupon were reused, not solved again.

Each model comparison uses its own fine frequency band. Do not compare the two cases by an intersection-only score: the physical arm change moves the resonances to different bands. Sorted dip minima are not field-identified modes and their splitting is not automatically 2g. Q3D shunt subtraction, Fast-sweep interpolation and unmodelled tee/bridge effects remain limitations.

[Spectra](physical_arms_s21.png) · [Exact results and raw-run sources](summary.json) · [Run provenance](run.json)

The parent controller was started with an earlier snapshot that attempted a common-fine-window ranking after the solves. That ranking, if present, is not used here because its overlap can exclude the resonances. This report supersedes that postprocessing; the actual HFSS raw data are unchanged.
