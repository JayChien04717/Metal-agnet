# Resonator simulations

Run from the repository root in qmetal with AEDT 2021 R2. Both commands default to the opposing double-U geometry.

```powershell
conda run -n qmetal python -B Design7/simulation/resonator/readout_ab.py q3d
conda run -n qmetal python -B Design7/simulation/resonator/readout_ab.py hfss
```

The runner records the rendered geometry, compact ground/substrate domain, boundary treatment, solver settings, raw exports and convergence in independent run folders. Add `--dry-run` for a geometry and ground-strap preflight. Analyze an explicit HFSS run with `readout_fit.py <hfss_run>` in the same environment. Compare it with the unfitted Q3D-based circuit using `readout_circuit.py <q3d_run> <hfss_run>`.

## Approved v4 layout

Pass the [approved v4 geometry input](../../results/resonator/20261118_140229_retune_v3_1daa0c12/raw/geometry_retune.json) to `readout_ab.py` via `--geometry`. Run Q3D and inspect the Maxwell nets, adaptive error and filter–resonator mutual first. Running `readout_circuit.py <q3d_run> --q3d-pass 15` freezes the π-circuit inputs and predicts J and both poles before HFSS. The model uses pass-specific Q3D mutuals and interpolates the 60/80 µm hammer coupon shunts at the approved head span.

Then run `readout_ab.py hfss --geometry <approved_geometry> --upstream <prediction_run> --wide 6.55 7.25 --fine 6.70 7.08`. Use `readout_circuit.py <q3d_run> <hfss_run> --frozen <prediction_run>` to compare the complex spectra with the recorded circuit values. The comparison applies the known 50 Ω port-polarity convention.

## Resonator-only eigenmodes

Run `readout_eigen_audit.py <v3_geometry> <v3_driven_run> <v3_prediction_run>` to reproduce the six v3 geometry comparisons. Add `--case isolated --min-freq 3.0` for the lower-frequency search on the original resonator. `readout_eigen_compare.py <eigen_run>` reports the cases and calculates the differences against their ideal circuit references. Add `--fields` to export stored gap fields for the no-tee case.

## Aligned short-leg draft

In qmetal, run `readout_retune.py <v3_sizing_run> <v3_frozen_circuit_run> <v3_s21_fit_run>` with the three explicit upstream directories from the [retune run record](../../results/resonator/20261118_140229_retune_v3_1daa0c12/run.json). It sizes a draft using circuit sensitivities and builds the layout without launching AEDT. Pass its `raw/geometry_retune.json` to `Design7/layout/design7_layout.py --figures --geometry` to reproduce the previews.
