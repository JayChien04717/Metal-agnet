# Design4 hammer: fixed-length Q3D / HFSS comparison

Exploratory simulations completed on 2026-11-16. These runs are not an adopted fabrication version. No filter or resonator length, phase velocity, capacitance, phase offset, or amplitude was fitted. The user authorized hammer gap/size and CPW ground configuration changes while keeping line lengths fixed.

## Findings

There is no demonstrated single best matching geometry in this four-candidate comparison. The dense-ground candidate minimizes the integrated complex S21 error, but has the largest mismatch between the two dip frequencies. A 100 um hammer height gives the smallest dip-frequency and splitting errors, but a worse integrated complex error than the baseline. None has yet demonstrated a converged, quantitatively matching coupled spectrum.

All scores below use the same 5.793–6.179 GHz interval and 25 kHz grid, with the deterministic port-polarity correction described below. RMS is an absolute dimensionless S-parameter difference, not a percentage. Dip pairs are sorted spectral minima; mode localization was not measured, and the splitting must not be interpreted as 2g at unknown detuning.

| Gap / head height (um) | Ground configuration | Complex S21 RMS | Dip-frequency RMS (MHz) | Splitting error (MHz) |
|---|---|---:|---:|---:|
| 42 / 65, baseline | auto bonds + hammer/tee straps (`all`) | 0.031969 | 31.682 | 37.400 |
| 42 / 65 | explicit CPW straps, nominal 200 um pitch (`cpw`) | **0.030480** | 52.940 | 66.550 |
| 20 / 65 | `all` | 0.042587 | 40.522 | 44.400 |
| 42 / 100 | `all` | 0.039701 | **22.915** | **24.550** |

For the 42 / 100 um candidate, HFSS dips are 5.874950 and 5.928550 GHz, versus circuit predictions 5.906575 and 5.935625 GHz. Their splitting is 53.600 versus 29.050 MHz. Thus even the best dip-frequency candidate is not a successful spectrum match. The narrow and broad features also occur in different orders in the model and HFSS; matching sorted minima alone does not establish mode correspondence.

Corrected [ranking and exact sources](../results/resonator/20261116_072453_hammer_verified_conventions_1c267037/parameter_report.json), [spectra](../results/resonator/20261116_072453_hammer_verified_conventions_1c267037/parameter_report.json), and [solver/port audit](../results/resonator/20261116_072453_hammer_verified_conventions_1c267037/parameter_report.json).

## Geometry and model

The experiment explicitly uses the hammer version from `layout/purcell_6500.py`, not the current generic `PurcellResonator` component, which has no hammer. Fixed filter/resonator lengths are 4.535956670430899 / 4.852581441726434 mm. CPW width/gap is 16 / 9 um, line spacing 400 um, substrate thickness 500 um, epsilon_r 11.5. Hammer head width is 25 um, ground clearance 16 um, corner radius 2 um. The full channel includes both Cin and Cp. Simulation bond wires are PEC; this is not a fabricated bridge design.

The circuit uses the same physical lengths and enclosing-capacitor coupon in all cases. Its quasi-static line model uses epsilon_eff=(1+epsilon_r)/2, Z0=49.89486 ohm, velocity=119916983.2 m/s. HFSS ports remain 50 ohm. Cin/Cp mutual capacitance is 39.5307 fF, with frame/plate ground terms 118.03591 / 17.10815 fF. These local Cin/Cp endpoint contributions have not been deembedded.

The hammer is a pi network. Mutual capacitance is minus the Maxwell off-diagonal; ground terms are row sums of the signal-only submatrix, excluding the reference ground. The extracted deembedded values (mutual / ground A / ground B, fF) are:

- 42 / 65 um: 1.45200 / 2.392485 / 2.410975.
- 20 / 65 um: 2.82917 / 3.361235 / 3.349655.
- 42 / 100 um: 2.52294 / 6.056165 / 6.066555.

The dense-ground case reuses the same local hammer coupon. The local electrostatic coupon does not model the full channel's finite bridge inductance or bridge loading; this is an explicit model limitation, not an exact equivalence of full Q3D/HFSS geometries.

## CPW extension and reference planes

The supplied excerpt says to extend CPW sufficiently far, but gives no numerical minimum. We tested 50, 100, 200, 400, and 800 um leads per side in a fixed 5.2 x 4.3 mm domain. The background reference is a continuous straight CPW with identical outer endpoints. Half its ground capacitance is subtracted from each hammer signal's ground term. The remaining pi network is placed at the hammer centre and the physical tee-to-centre arms are represented as transmission lines.

The longest two cases were continued on their existing Q3D adaptive setups to a maximum of 20 passes / 0.1% requested error. Their original eight-pass histories were preserved. Final values:

| Lead per side (um) | Mutual (fF) | Excess ground A / B (fF) |
|---|---:|---:|
| 400 | 1.43584 | 1.634235 / 1.643775 |
| 800 | 1.45200 | 2.392485 / 2.410975 |

Mutual changes by about 1.1%; the excess shunts still differ by about 0.75 fF. At 800 um the last-pass subtraction sensitivity is about 0.56 / 0.58 fF. This is a numerical sensitivity indicator, not a rigorous error bar. **800 um is the longest tested choice, not a demonstrated converged extension length.**

For gap20, the lead is 811 um to retain the same outer endpoint at +/-846 um and reuse the exactly matching reference. Geometry comparisons verified this. Only new hammer coupons were solved for the changed gap and height.

Original sweep: [run](../results/resonator/20261115_234035_hammer_cpw_deembed_eff4e52d/parameter_report.json). Continued meshes and final CPW data: [run](../results/resonator/20261116_000326_hammer_cpw_refine_2d2775df/parameter_report.json).

## What deembedding changed

For the baseline's own fine-sweep interval, after correcting port polarity:

| Circuit treatment | Complex S21 RMS |
|---|---:|
| Diagnostic omitting hammer's CPW arms | 0.077889 |
| Raw 10 um coupon with correct physical pin planes | 0.039235 |
| Deembedded 800 um coupon with centre reference planes | 0.031336 |

The raw coupon uses tee-to-pin arms of 0.144 mm each; the deembedded model uses tee-to-centre arms of 0.200 mm each. This retains the same physical channel and avoids counting CPW capacitance twice. Explicitly representing the arms moves predicted dips by roughly 130 MHz relative to the omitted-arm diagnostic. Deembedding then improves complex RMS by approximately 20%, but it does not align both dips or establish that truncation is the sole cause of the mismatch.

Corrected [baseline analysis](../results/resonator/20261116_072510_hammer_model_audit_1c44aec2/parameter_report.json) and [figure](../results/resonator/20261116_072510_hammer_model_audit_1c44aec2/parameter_report.json). This analysis reused raw EM exports; it did not rerun solvers.

## Convention corrections and convergence limits

Qiskit Metal `create_ports` calls pyEPR `Rect.make_lumped_port(axis)`, which places both integration lines in the positive global coordinate direction. The input/output pin normals are (-1,0) and (+1,0). Mapping the nodal circuit's signal-ground convention to these ports multiplies S21 by -1, with no frequency-dependent phase or fitted constant. Renderer source excerpts, rebuilt geometry checks and port normals are saved in the corrected ranking run. All magnitude spectra are unchanged.

The original baseline/search complex-error summaries omitted this polarity conversion. Their raw solver exports remain usable, but their complex-error rankings are superseded by the corrected runs linked above. The early `20261116_005237_hammer_model_audit_3bea5c33` analysis also predates this correction.

The original CPW sweep's `caps.json` last-pass relative-change calculation used a denominator floor inappropriate for matrices exported in F. The corrected analyses recompute changes from the original all-pass CSVs. Final capacitance values in fF were unaffected. These different audit metrics must not be confused with the solver's own convergence criterion.

HFSS retained its original three-pass, 10 um seed, first-order setup. Last Delta S values are 0.046356, 0.044907, 0.033703 and 0.036799 for baseline, dense-ground, gap20 and height100 respectively, all above the requested 0.01. Fine-sweep power-balance errors are at most 5.74e-4; this does not establish mesh or interpolation convergence. No direct-frequency solves or field-based mode assignments were performed.

Q3D and HFSS used shared geometry definitions, but not a shared solver mesh. Cross-solver mesh reuse was not demonstrated. Each HFSS wide/fine pair uses the identical adaptive Setup and unchanged adaptive-history record. Changed geometry requires its own HFSS mesh. Identical reference and Cin/Cp Q3D results were reused after geometry checks. The historical HFSS design sought for reuse was absent, so all four full-channel spectra here were newly simulated.

Remaining model limitations include distributed electric/magnetic coupling, tee and bend discontinuities, finite ground-bridge loading, Cin/Cp endpoint deembedding, and the quasi-static line propagation model. Numerical convergence and these model omissions prevent a validated optimum claim. Increasing head height is the best tested direction for the two-dip locations; denser ground is the best tested integrated-error case. Neither should be released as a confirmed fit.

## Reproduction

Use the Design4 `qmetal` environment. Entry points are in [simulation/resonator](../simulation/resonator/):

1. `hammer_deembed.py`: extension sweep and raw/reference Q3D extraction.
2. `hammer_refine.py`: continue the explicitly selected parent setups.
3. `hammer_fixed_s21.py`: fixed-length HFSS and circuit comparison.
4. `hammer_search.py BASELINE_RUN REFERENCE_RUN`: bounded geometry search; launches sequential EM solves.
5. `hammer_reanalyse.py BASELINE_RUN REFERENCE_RUN`: offline baseline diagnostics.
6. `hammer_finalize.py SEARCH_RUN`: offline polarity correction for the explicitly recorded pre-correction search. It rejects already corrected data to prevent double flipping.

Every run records inputs, code snapshots, raw exports, environment and upstream references. The search is [20261116_004926_hammer_geometry_search_32bc0e37](../results/resonator/20261116_004926_hammer_geometry_search_32bc0e37/parameter_report.json). Its progress record names the four individual HFSS runs and two candidate Q3D runs. No GDS was generated and no final-selection or release-confirmation state was changed by this experiment.

## Follow-up: actual 400 / 800 um arms in HFSS

The user subsequently requested actual full-channel 400 / 800 um CPW arms and confirmed that the two main resonator lengths remain fixed. Two new HFSS geometries were solved, using head-outer-edge to main-line tap length per side as the definition. This differs from changing only the local Q3D extraction boundary. The straight-arm construction changes main-line centre separation to 0.892 / 1.692 mm, so distributed inter-line coupling can also change.

| Actual arm per side | HFSS dips (GHz) | Deembedded circuit dips (GHz) | Last HFSS Delta S |
|---|---|---|---|
| 400 um | 5.662875 / 5.727450 | 5.724725 / 5.740100 | 0.019345 |
| 800 um | 5.323525 / 5.420950 | 5.394425 / 5.422025 | 0.022242 |

Both cases retain gap42 / height65 / width25 um, CPW16/9 um, main lengths4.535956670430899 / 4.852581441726434 mm, and the existing ground policy. Both remain above the Delta S target0.01. Fine-sweep adaptive histories are unchanged; wide/fine mesh reuse is verified within each case. Maximum two-port power residuals are9.25e-6 and1.12e-3, respectively. No direct-frequency or field-localization validation was performed. These are exploratory EM spectra, not confirmed mode assignments or converged final values.

At400 um, local-band complex RMS is0.027094 for the raw full-arm coupon and0.032024 after deembedding. At800 um, it improves from0.030548 to0.015552. The800 um high-frequency dip agrees within1.075 MHz, but the low-frequency dip differs by70.900 MHz. Thus the physical extension does not by itself resolve the two-mode mismatch. These RMS values use each case's own fine band, not a common-frequency ranking.

[Completed report and definitions](../results/resonator/20261116_082546_hammer_physical_arms_report_ea38aecf/README.md) · [actual spectra](../results/resonator/20261116_082546_hammer_physical_arms_report_ea38aecf/physical_arms_s21.png) · [exact results and sources](../results/resonator/20261116_082546_hammer_physical_arms_report_ea38aecf/summary.json).

Reproduce with `hammer_physical_arms.py REFERENCE_RUN BASELINE_RUN`, then `hammer_arms_report.py PARENT_RUN`, in the Design4 qmetal environment. The parent is `20261116_074123_hammer_physical_arms_fe38e68b`; it references the existing refined Q3D run and baseline Cin/Cp extraction. Only the HFSS geometries were newly solved. Its earlier intersection-only postprocessing is superseded by the report above. The initial report attempt `20261116_082508_hammer_physical_arms_report_822c6253` failed during summary finalization and is not the report to use; raw EM data were unaffected.

## Follow-up: each arm 50 um, head height 30 um

The latest user instruction supersedes the proposed 100 um arms: each actual arm is 50 um, measured from the head outer edge to the main-line tap centre. Head height is 30 um, width remains 25 um, and the pad gap must be below the original 42 um. The target is to preserve Purcell–readout coupling, not the old local mutual capacitance. Main lengths remain 4.535956670430899 / 4.852581441726434 mm. The original gap42/height65/arm154 design is the comparison reference.

The short-arm construction reduces main-line centre spacing; it overrides the older 400 um spacing preference under the explicit user geometry request. Meander pitch is unchanged. Two local hammer ground wires use 12 um diameter: their 17 um landing interval cannot fit the previous 15 um wire with the existing edge margins. Other ground definitions remain unchanged. Geometry checks passed, but this is not fabrication signoff.

The original HFSS setup was continued without rebuilding geometry. Its first three adaptive rows were retained and the eighth pass reached Delta S 0.0087498, below the 0.01 target. [Refined baseline](../results/resonator/20261116_105611_hammer_hfss_refine_4229713c/parameter_report.json) and [shared-pole analysis](../results/resonator/20261116_112957_hammer_shared_poles_84ccce84/parameter_report.json) give conditional effective J = 10.576423 MHz. The older three-pass estimate was 10.199343 MHz; old dip locations and circuit-error rankings are historical numerical results, not converged baseline values.

The arm50/head30/gap20 candidate was actually solved with Q3D and HFSS. Local mutual C is 1.071230 fF; subtracting half of the matched straight-CPW reference gives signed excess ground values -0.355865 / -0.355875 fF. These are reference differences, not isolated physical negative capacitors. HFSS continuation retained its first three mesh-history rows, but pass8 Delta S = 0.3125 still fails 0.01. Conditional J changed from 17.673333 MHz at pass3 to 18.136306 MHz at pass8, approximately 71.5% above the refined original reference. This candidate does not demonstrate preserved coupling. [Geometry, actual spectra and comparison](../results/resonator/20261116_113029_hammer_short_report_0a733609/parameter_report.json).

Couplings above come from common complex poles fitted to S11/S21, followed by inversion of H/h = [[fp+i*kappa_p/2, J], [J, fr]], assuming real reciprocal J and negligible intrinsic readout loss. They are not half the visible dip separation or independent HFSS scalar outputs. A synthetic known-J recovery check passed; polynomial-background degree sensitivity is small but does not bound EM error. No field-mode assignment or direct-frequency validation was performed.

Reusing frozen Q3D circuit predictions on identical frequency grids, refined-baseline complex S21 RMS is 0.058116 raw / 0.050372 deembedded ([analysis](../results/resonator/20261116_113142_hammer_refined_model_063b6d5f/parameter_report.json)). Gap20 gives 0.051916 raw / 0.049208 deembedded ([analysis](../results/resonator/20261116_105856_hammer_refined_model_e1f5df8e/parameter_report.json)). Different fine bands prevent treating these as a common-band ranking. No model parameters or line lengths were fitted. Analysis `20261116_113108_hammer_refined_model_0193084d` is explicitly invalidated because it referenced historical circuit data without the port-polarity correction; its raw files are retained and the corrected analysis is the baseline link above.

Reproduction: `hammer_short_candidate.py BASELINE_RUN --gap GAP --max-passes 8`; `hammer_hfss_refine.py FULL_RUN --max-passes 8`; `hammer_pole_fit.py FULL_OR_REFINED_RUN`; `hammer_short_report.py CANDIDATE_RUN BASELINE_POLE_RUN CANDIDATE_POLE_RUN`. The eight-pass budget is an explicit validation override; the normal HFSS runner still defaults to three passes. Q3D and HFSS use shared geometric definitions, not a demonstrated shared solver mesh.

## Current search constraint: fixed 190 um main-line centres

The user subsequently confirmed fixed main-line centre spacing190 um and head height30 um. Head gap, head width and ground configuration may vary; each symmetric arm is derived as (190 - gap - 2*head width)/2 um. This supersedes the fixed50 um arm constraint and the earlier requirement to keep gap below42 um. Both main physical line lengths remain fixed. The aim is to retain the original effective coupling and independently improve the Q3D circuit versus HFSS spectrum; a pole fit to HFSS alone is not evidence that the Q3D model agrees.

The width25/gap40/arm50 case has actual Q3D mutual C0.488050 fF and an eight-pass conditional J12.667352 MHz, versus refined original10.576423 MHz (+19.77%). Its Delta S0.12097 fails0.01; continued mesh validation is recorded separately. [Eight-pass report](../results/resonator/20261116_125030_hammer_short_report_4a759351/parameter_report.json). No final parameter selection follows from this result.

One prepared candidate is width20/gap51/arm49.5 um, centre spacing190 um, height30 um, existing `all` ground configuration and12 um local bridge diameter. [Geometry-only check](../results/resonator/20261116_131650_h30_g51_w20_a49.5_27f0d106/parameter_report.json) passed; this dry run contains no EM result. Width20/gap52/arm49 reached the strict bridge-landing boundary and failed the geometry check, so it was not simulated. A preceding long-named width12 dry-run attempt failed during Windows path-length-limited input snapshotting, before geometry or EM execution; its manifest records preparation_failed. Shorter run tags resolve that preparation issue without changing shared run-recording rules.

The subsequent [fixed-circuit pole diagnostic](../results/resonator/20261116_145728_hammer_model_poles_24a3d66f/parameter_report.json), on the gap40/width25 eight-pass spectrum over5.986–6.25 GHz, gives raw/deembedded model J2.871280/2.883346 MHz versus conditional HFSS J12.667352 MHz. Deembedded bare-frequency errors are approximately -29.78/-21.47 MHz. These are unchanged circuit predictions, not fitted capacitances or line lengths. The discrepancy means tuning HFSS coupling alone cannot establish model agreement. Ground-return/distributed effects are hypotheses, not demonstrated causes. The reused Cin/Cp coupon also reports nonconvergence and needs separate continuation before its uncertainty is excluded.

Consequently, the next diagnostic priority is unchanged gap40/width25/height30/arm50 geometry with a controlled CPW-ground configuration change, after completing the active HFSS refinement and refining Cin/Cp. The width20/gap51 candidate remains deferred. Local Q3D coupons do not include the full channel's three-dimensional ground bridges; reusing them for this diagnostic does not establish bridge-inclusive electrostatic equivalence or a common Q3D/HFSS mesh. The existing `cpw` option replaces automatic bonds with explicit CPW-segment bridges plus the hammer/T-junction bridges; this is a ground-configuration comparison, not a single-bridge perturbation.

The [CPW-ground diagnostic layout](../results/resonator/20261116_150518_h30_g40_w25_a50_43c7bfb3/parameter_report.json) passed geometry and bridge-foot/etch checks with fixed dimensions; this is a dry run only. An independent [original-layout model pole comparison](../results/resonator/20261116_151024_hammer_model_poles_d3df9815/parameter_report.json) gives raw/deembedded J6.571264/8.105286 MHz versus refined HFSS10.576423 MHz. Its deembedded Purcell bare-frequency error is -29.94 MHz, similar to the short-arm case, while the readout error is +1.02 MHz. Thus part of the discrepancy predates the short-arm layout, and the larger mismatch at190 um spacing must be investigated separately. Pole parameters remain conditional on the two-mode model; the two analysis bands differ, so their whole-band RMS errors are not a common-band ranking.

The gap40/width25 continuation has now completed: [eleven-pass HFSS](../results/resonator/20261116_125008_hammer_hfss_refine_f5653eec/parameter_report.json) retains its earlier mesh history and reaches Delta S0.0052015 at1,686,549 tetrahedra. Its [updated comparison](../results/resonator/20261116_152648_hammer_short_report_8befc087/parameter_report.json) gives conditional J12.752274 MHz, +20.57% versus the original reference and +0.67% versus pass8. The [unchanged-circuit analysis](../results/resonator/20261116_152545_hammer_model_poles_abeb6ef8/parameter_report.json) still predicts deembedded J2.883346 MHz, with bare-frequency errors -30.90/-23.62 MHz. The substantial mismatch survives this mesh continuation.

A separate [spectrum audit](../results/resonator/20261116_152545_hammer_spectrum_audit_d93bb40a/parameter_report.json) finds maximum two-port power residual5.26e-8. The archived final-adaptive solution at6.08 GHz agrees with the Fast sweep to absolute complex S21 error1.30e-8 (S11 error2.27e-9). This validates one directly solved frequency on the final mesh; dedicated solves at both resonance minima and field-mode identification remain unperformed. The report's statement about absent direct-frequency validation should be read with this separate one-point audit. Fast sweeping was memory/out-of-core limited and completed before any subsequent Ansys solve was started.

At the user's request, plots now explicitly identify the solver/model. [Driven versus Q3D circuit](../results/resonator/20261116_152952_hammer_refined_model_7e4361c7/parameter_report.json): black solid = HFSS Driven; blue dashed = raw Q3D capacitance circuit; orange dash-dot = deembedded Q3D capacitance circuit. Q3D does not directly produce these S21 curves: the extracted capacitances feed the fixed transmission-line circuit. The lower panel is absolute complex S21 difference from Driven. [Original versus short-arm geometry comparison](../results/resonator/20261116_153038_hammer_short_report_88ef9724/parameter_report.json) labels both curves as HFSS Driven. These replots reuse existing EM data and do not introduce new solver results.

The figure above freezes the earlier Cin/Cp inputs. A [20-pass Cin/Cp continuation](../results/resonator/20261116_152626_hammer_cp_refine_cbfebf1b/parameter_report.json) preserves its first8 passes and changes mutual C39.53070 to41.25192 fF, with ground values122.27552/17.50270 fF. Native Delta0.40053% still fails the0.1% target despite the last matrix change being approximately0.142%. The [circuit recalculation using these inputs](../results/resonator/20261116_153202_hammer_model_poles_48760dfb/parameter_report.json) gives deembedded J2.885092 MHz versus HFSS12.752274 MHz. Thus this observed capacitance correction does not explain the coupling discrepancy; it increases the Purcell frequency error. This is a sensitivity comparison with still-unconverged Cin/Cp inputs, not their final validation.
