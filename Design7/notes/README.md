# Design7 layout notes

## Experiment objective

The 2026-11-18 user request defines Design7 as a test of layout-driven agreement between the standard Q3D-to-circuit workflow and driven HFSS S21. The Design6 side-by-side CPWs introduced distributed mutual capacitance, direct feed–resonator coupling and ground-return sensitivity. Design7 tests whether separating the long CPW runs while keeping a compact short hammer makes those contributions small enough for the ordinary circuit model to predict S21 without a Design7-specific capacitance formula, fitted capacitance multiplier or added parasitic compensation. [Design6 methods and evidence](../../skill/purcell_resonator_design.md)

Freeze the as-drawn dimensions, CPW parameters and Q3D-derived capacitance inputs before comparing with HFSS. Check both resonance frequencies, linewidths and the complex S21 trace over the readout band. Compare J only when the poles and their linewidths support the inversion. Record actual adaptive error and sweep resolution alongside discrepancies. No numerical agreement tolerance was specified. A full-geometry Maxwell mutual includes every electrostatic path, so its value alone cannot establish that distributed coupling is negligible.

## Fixture and simulation scope

The user sketches place the enclosing feed coupler above the filter, with opposing hammer electrodes between two vertical stems. The compact revision folds the filter shorted tail to the right and the resonator shorted tail to the left. The downward resonator stem leads to the readout open end. The planar schematic was interpreted as a single-channel fixture with 50-ohm feedline port locations. The branch tees and folds are draft geometry.

The baseline takes Q1 line lengths and hammer electrodes from [Design6 r10 wafer_candidate](../../Design6/design.json). The loaded Q1 resonator length is 4.053848 mm. The no-claw fixture adds [0.100 mm CLAW_DL](../../Design6/layout/design6_layout.py) to both layouts. The 1.20 mm and 1.30 mm tap distances, the enclosing feed capacitor, and the 16/9 µm CPW come from [the existing channel builder](../../Design4/layout/purcell_6500.py). The tee-arm dimensions and routing in the opposing variant are Design7 draft inputs. Consult [geometry.json](../parameters/geometry.json) for dimensions.

The user requested Q3D and HFSS S21 for the double-U geometry on 2026-11-18 and directed that existing data serve as the reference. Q3D opens the feed and shorted-tail ends to identify the three signal nets. The Maxwell filter–resonator mutual contains the hammer and all other electrostatic paths of the drawn geometry. Driven HFSS uses the original shorted ends, 50 Ω lumped feed ports and PEC ground straps. Check tap loading, bare detuning, linewidth resolution, and adaptive convergence before interpreting a two-mode estimate of J. A near-resonance retune would constitute a separate geometry revision.

The comparison uses the compact 6.2 × 4.4 mm ground and silicon domain centred at (0, −1.55) mm. The feed port planes are at x = ±2.5 mm, with approximately 0.6 mm perimeter clearance. Both solvers retain the as-drawn metal and a 0.5 mm silicon thickness. The smaller domain is a fixture boundary change. Compare with historical results only after checking the applicable reference planes and ground extent.

The unfitted prediction follows the existing full-Q3D and local-coupon π convention. It uses the new fixture's Q3D feed–filter and filter–resonator Maxwell mutuals for C_p and C_j. Nominal shunts come from the recorded Design6 coupons for matching electrode geometry. The hammer shunts retain the standard straight-line deembedding. The circuit uses the drawn branch lengths and conformal 16/9 µm CPW values. Direct feed–resonator mutual is reported separately and omitted from this test of the two-capacitor circuit. [Design5 model convention](../../Design5/simulation/resonator/theory5.py) · [Design6 nominal coupon record](../../Design6/results/resonator/20261116_221901_readout_caps_0e896b6c/summary_out.json)

## Execution history

The first Q3D attempt used an interior route pin as the right feed opening and combined the feed with ground. Its raw matrix failed the net and row-sum checks. The corrected solver runner now checks both port locations and outward normals against the drawn endpoint components.

The corrected compact-ground HFSS job ran ahead of the valid Q3D extraction. Its coupling had therefore not been checked when the sweep started. The comparison uses the subsequent valid compact-ground Q3D data as independent model inputs.

For the same 10 µm signal seed, the recorded HFSS initial meshes had about 249,100 tetrahedra with the 10 × 10 mm ground and 254,845 with the 6.2 × 4.4 mm ground. The compact case reduces the ground area by 73%, with no initial mesh-count saving in this pair of runs.

The active Q3D job had passed 15 when the user capped the study at 15 passes. AEDT's Clean Stop completed the current pass 18 and preserved adaptive data. Passes 14 and 15 were exported without another solve. The comparison selects pass 15, and the solver configuration now caps subsequent runs at 15. [Pass-15 export](../results/resonator/20261118_110609_q3d_pass_export_2953be42/run.json)

## Outcome

The extracted direct feed–resonator mutual is 0.01317 fF. The 0.90397 fF filter–resonator mutual includes hammer and distributed paths. A comparison of the same signal geometry and ground domain gives a 33.75 MHz broad dip offset and a 3.62% broad linewidth difference between the unfitted Q3D-based circuit and HFSS. The phase overlay converts the predetermined sign of one HFSS port integration line to the circuit convention. No phase delay, capacitance, line length or multiplier was fitted. The Design6 nominal shunt coupon also had incomplete convergence.

The two-capacitor model reproduces the broad linewidth to within 0.68 MHz in this comparison. The spectrum still has a frequency offset and the 0.25 MHz fine fast sweep resolves only the broad pole. The current evidence does not establish agreement for the resonator mode or a reliable dynamic J. [Model comparison summary](../results/resonator/20261118_110804_q3d_circuit_comparison_6f5e972d/summary.json)

A separate circuit-only reanalysis estimates J/2π = 6.15 MHz from the two bare driving-point susceptance slopes at the hammer. The modeled bare frequencies are 5.840 and 6.910 GHz, a 1.069 GHz detuning. The resulting resonator-like mode has an external linewidth of about 1.17 kHz in the model with zero intrinsic loss. This estimate uses the full Q3D mutual localized at the hammer and excludes distributed magnetic coupling. Design6 r10 loaded HFSS fits gave 15.96 and 15.25 MHz. [Design7 circuit estimate](../results/resonator/20261118_111507_q3d_exchange_estimate_652407ec/summary.json) · [Design6 loaded HFSS](../../Design6/README.md#loaded-readout-s21)

## 15 MHz draft

Revision `d7_readout_15mhz_v3` increases the transverse span of each hammer head from 37.5 to 74.2 µm. The 25 µm axial thickness, 40 µm inter-head gap and hammer axis preserve the compact end-to-end coupling arrangement. Shortening the outer return leg of the right-side filter U by 765 µm moves its grounded endpoint from approximately y = −0.477 mm to −1.242 mm and reduces the filter's total CPW length to 3.1154 mm. The feedline, resonator, tees, and fold coordinates use the v2 dimensions.

The approximate sizing interpolated the recorded 60 and 80 µm Design4 local hammer results, including the deembedded shunts, and carried forward the v2 full-minus-local mutual of 0.240 fF as a first estimate. The resulting sizing input was 1.915 fF. The circuit calculation targeted 15 MHz exchange with the bare modes together near 6.873 GHz. [Sizing assumptions and outputs](../results/resonator/20261118_111951_size_15mhz_candidate_5e511d71/run.json)

The user approved this layout for EM testing. The compact-ground Q3D run completed 15 passes before HFSS started. Its filter–resonator mutual is 1.96890 fF, direct feed–resonator mutual is 0.01400 fF, and feed–filter mutual is 41.73450 fF. The filter–resonator mutual changed by 0.1351% from pass 14 to 15, with global ΔC of 0.41718% against a 0.1% target. Before HFSS, the calculated mutuals, dimensions and coupon shunts were frozen in the prediction run. The resulting model predicts J/2π = 15.42 MHz and a −0.04 MHz bare detuning. [Q3D run](../results/resonator/20261118_112242_q3d_ab_ffee939e/run.json) · [Frozen circuit inputs](../results/resonator/20261118_112740_q3d_circuit_prediction_6ba5a9a7/run.json)

Three HFSS adaptive passes gave ΔS = 0.028285 against a 0.01 target. The 0.25 MHz fine fast sweep resolves two poles. Their 0.647 and 28.072 MHz fitted linewidths give J/2π = 14.154 MHz and κ_p/2π = 28.718 MHz under the one-lossy-mode inversion. The 1 MHz wide fast sweep gives frequencies and inversion consistent with the fine sweep, with the narrow linewidth below the wide grid-spacing criterion. The fine result is 5.64% below the 15 MHz target. It remains conditional on the two-mode model and a 90.08 MHz detuning, which exceeds twice the fitted exchange.

The frozen model exchange exceeds the HFSS fit by 1.27 MHz. The same two-pole fit on the frozen model gives bare f_p/f_r = 6.87271/6.87238 GHz, versus 6.90753/6.81745 GHz in HFSS. In particular, the HFSS filter is 34.82 MHz higher and the HFSS resonator 54.93 MHz lower than modeled. The spectrum therefore has a different splitting and linewidth distribution. The direct complex RMS difference is 0.353 on the 6.70–7.08 GHz grid using the predetermined port polarity. The comparison used a single geometry and domain, confirmed by matching geometry snapshot hashes, and retained all pre-HFSS circuit inputs. [HFSS run](../results/resonator/20261118_112807_hfss_ab_209bf4a4/run.json) · [Complex overlay](../results/resonator/20261118_114138_q3d_circuit_comparison_b526a6b6/s21_q3d_vs_hfss.png)

### Resonator frequency audit

Independent fits of S11, S21, S12 and S22 all give a bare resonator at 6.81745 GHz on the same HFSS mesh. The narrow hybrid pole is at 6.81533 GHz, only 2.12 MHz below the inferred bare frequency. The exact bare circuit root is 6.87249 GHz. Both EM inputs have the same geometry hash and the model's open-to-short CPW path equals the layout length of 4.153848 mm. The saved AEDT sweep path for the v3 `trace_res_tail` includes all three 99 µm fillets.

In the frozen circuit, the 55.04 MHz difference is equivalent to adding 34.49 µm to the resonator's shorted arm or 7.41 fF at the resonator hammer plane. These are independent sensitivities. A uniform 0.805% velocity reduction would match the resonator and would move the model filter down to 6.8210 GHz, farther from its HFSS value of 6.90753 GHz. Also, the full Q3D resonator's total static ground capacitance is 18.11 fF below the transmission-line-plus-coupon total. The whole-net ground capacitance alone cannot locate local loading or determine modal phase. The precise contribution of the tees, short, open, bends, and hammer end remains undetermined at the three-pass HFSS resolution. [Audit](../results/resonator/20261118_120541_bare_frequency_audit_079373e5/summary.json)

### Resonator-only HFSS eigenmodes

Six v3 cases retain the 4.153848 mm open-to-short resonator path, the 6.2 × 4.4 mm domain, and common resonator-side straps. The first search used a 5.5 GHz floor, two modes, a 10 µm metal seed, and six adaptive passes. The companion modes are above 21 GHz. Each comparison uses its own ideal transmission-line reference so that removing the open stub does not get counted as an electromagnetic tee residual.

| Change | HFSS frequency change (MHz) | Ideal circuit change (MHz) | Excess HFSS change (MHz) |
|---|---:|---:|---:|
| Ground the opposing hammer electrode instead of omitting it | −0.63 | −14.80 | +14.17 |
| Replace the plain CPW stub with the hammer and grounded mate | −21.38 | −36.60 | +15.23 |
| Add the tee and its plain open stub | −304.59 | −308.13 | +3.54 |
| Increase the resonator open-end gap from 9 to 45 µm | +2.53 | 0 | +2.53 |
| Flare the last 75 µm of the shorted route to a 28 µm terminal width | +28.39 | 0 | +28.39 |

The original isolated resonator is at 6.88099 GHz, with its grounded hammer mate at 6.88036 GHz. The latter is 7.87 MHz above the corresponding frozen circuit frequency and 62.91 MHz above the bare frequency inferred from the complete driven assembly. The last isolated-mode adaptive steps are 3.14–6.84 MHz across the six cases. The original isolated frequency also increased by 25.64 MHz from adaptive pass three to six. This identifies local short-contact sensitivity and a separate discrepancy between the resonator-only and complete driven setups; the measurements do not separate the larger assembly difference into environment and mesh contributions.

The 16-to-28 µm flare changes the short-region impedance over a finite length. Its 28.39 MHz shift quantifies that geometry change. The absolute frequency reduction of the original short remains undetermined. The open-gap shift is below the last-step scale. The supplemental gap-field phase fits depend strongly on fit window and hit the allowed open extension bound, so they were excluded from the endpoint attribution. Each original solver case exported its complete frequency-vs-pass table, then the wrapper failed while serializing a missing lossless linewidth. The comparison recovered the final frequencies from those tables without an EM rerun. The exporter now handles missing linewidths. [Recovered results and convergence plot](../results/resonator/20261118_132936_eigen_audit_fit_76ea3be0/eigenmode_comparison.png) · [Numerical audit](../results/resonator/20261118_132936_eigen_audit_fit_76ea3be0/summary.json)

The repeat search lowered the eigenmode floor from 5.5 to 3.0 GHz on the original isolated geometry, with an identical geometry snapshot hash and the same two-mode, six-pass budget. Its ordered modes are 6.88190 GHz and 21.53973 GHz. It found no lower mode between 3.0 and 6.88190 GHz. The resonator frequency differs by +0.912 MHz from the 5.5 GHz-floor result. The final maximum adaptive change was 0.10141% against the 0.1% target. The solver reports a lossless linewidth as unavailable. [Lower-floor search](../results/resonator/20261118_133004_resonator_eigen_audit_d61ef386/raw/isolated/result.json)

The layout reduces close side-by-side overlap. The clearance increase is geometric evidence.

## Aligned short-leg draft

The `d7_readout_aligned_v4` draft uses the [shared short-tail method](../../note/purcell_filter_method.md#quarter-wave-frequency-and-short-tail-trim) and retains the original sizing target of 6.872891 GHz. The local filter and resonator slopes are −1.658 and −1.609 MHz/µm. The draft adds the finite circuit frequency shifts to the assembled v3 HFSS bare frequencies of 6.907527 and 6.817451 GHz. This assumes local transfer of the v3 offsets from the three-pass driven fixture. The isolated eigenmode remains a separate diagnostic.

| Element | Local change | Draft dimension |
|---|---|---|
| Right filter U | Extend the outer straight return leg at the ground short by 22.31 µm | 3.137710 mm total CPW |
| Left resonator U | Shorten the outer straight return leg at the ground short by 32.96 µm | 4.120889 mm total CPW |
| Both hammer heads | Reduce transverse span by 1.93 µm per head | 72.22 µm span, 25 µm axial thickness, 40 µm gap |

The exchange estimate uses the v3 frozen circuit capacitances and the 60–80 µm coupon derivatives. It gives J/2π = 15.00 MHz. The v3 HFSS exchange of 14.15 MHz remains conditional on the strongly detuned fit. With the v3 offsets carried forward, the draft bare modes both land at 6.872891 GHz. The raw updated circuit roots are 6.837890 and 6.927926 GHz. Qiskit Metal checks found three connected signal nets, no inter-net short, 457.882 µm filter–resonator clearance outside the hammer region, and valid port planes and ground straps. The compact-ground footprint remains 6.2 × 4.4 mm. The user approved this layout for Q3D and HFSS on 2026-11-18. [Sizing assumptions and checks](../results/resonator/20261118_140229_retune_v3_1daa0c12/summary.json) · [Reproduction](../simulation/resonator/README.md#aligned-short-leg-draft)

### v4 compact-ground results

The driven HFSS filter moved down 35.26 MHz relative to v3 and landed 0.63 MHz below the common target. The resonator moved up 27.38 MHz and landed 28.06 MHz below the target. The absolute detuning improved from 90.08 to 27.43 MHz. The wide and fine fast sweeps give 13.74311 and 13.74315 MHz exchange on the same adaptive mesh. The fitted pole widths of 3.906 and 25.084 MHz both exceed twice the 0.25 MHz fine grid spacing. This check supports resolving the two linewidths on the exported interpolant. The inversion assumes one lossy mode, zero intrinsic readout loss, and real J. [v4 HFSS fit](../results/resonator/20261118_170017_opposing_s21_fit_2343c01c/summary.json)

The pass-15 Q3D prediction uses the static mutual at the hammer and gives J/2π = 14.97 MHz from the susceptance slopes. Applying the same two-pole inversion to the circuit spectrum gives 15.39 MHz versus 13.74 MHz in HFSS. The unadjusted Q3D model places the bare filter about 34.5 MHz below HFSS and the bare resonator about 83.1 MHz above HFSS. The direct amplitude and phase comparison retains the frozen circuit parameters and predetermined port polarity. The fine sweep used the established 6.70–7.08 GHz band and left the three-pass adaptive history unchanged. [Complex overlay](../results/resonator/20261118_170023_q3d_circuit_comparison_a4edb91d/s21_q3d_vs_hfss.png)
