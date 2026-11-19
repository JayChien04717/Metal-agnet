# Design7 Purcell filter and readout resonator

**Figures:** [fig/](fig/README.md) holds the whole chip ([`chip.png`](fig/chip.png)) and one close-up per part, rendered from the current GDS by `code/general/plot_design_figs.py`.

Design7 is a single-channel double-U simulation fixture. The user approved revision `d7_readout_aligned_v4` for Q3D and HFSS testing. The experiment tests whether a standard Q3D-derived circuit with inputs frozen before HFSS can predict complex S21 for the compact geometry. [Design record](design.json) · [Method and earlier results](notes/README.md)

## v4 layout

The filter and resonator U lengths are 3.13771 and 4.12089 mm. The resonator includes 0.100 mm of compensation for the omitted claw. The hammer heads have a 72.22 µm transverse span, 25 µm axial thickness and 40 µm separation. Both solvers use the 6.2 × 4.4 mm ground on 0.5 mm silicon. The sizing targets a common bare frequency of 6.872891 GHz and J/2π = 15 MHz. Layout (`design7_layout_opposing_geometry_retune.png`, removed 2026-11-19; in git history) · Hammer detail (`design7_layout_hammer_geometry_retune.png`, removed 2026-11-19; in git history) · [Approved geometry input](results/resonator/20261118_140229_retune_v3_1daa0c12/raw/geometry_retune.json) · [Trim assumptions](notes/README.md#aligned-short-leg-draft)

## Simulation results

Q3D completed 15 passes on the compact ground. The filter–resonator mutual is **1.90837 fF**, the feed–resonator mutual is 0.01407 fF, and the feed–filter mutual is 41.74969 fF. The filter–resonator entry changed by 0.149% from pass 14 to 15. The global adaptive ΔC was 0.38899% against the 0.1% target. [Maxwell matrix](results/resonator/20261118_163603_q3d_ab_b55f8759/raw/opposing_full/matrix_fF.csv)

| Quantity | Frozen Q3D circuit | HFSS S21 two-pole inversion |
|---|---:|---:|
| Bare filter (GHz) | 6.837805 | 6.872262 |
| Bare resonator (GHz) | 6.927951 | 6.844831 |
| Bare detuning f_r − f_p (MHz) | +90.15 | −27.43 |
| J/2π (MHz) | 14.97 | 13.74 |

The Q3D column uses the bare roots and the hammer-plane exchange calculation frozen before HFSS. The HFSS fit resolves both linewidths on the 0.25 MHz fine grid and gives a resonator 27.43 MHz below the filter. Its exchange is 1.26 MHz below the target. Three HFSS passes gave ΔS = 0.024178, acceptable under the shared driven criterion and above the configured 0.01 target. [Frozen prediction](results/resonator/20261118_164206_q3d_circuit_prediction_68f91ab7/summary.json) · [HFSS pole fit](results/resonator/20261118_170017_opposing_s21_fit_2343c01c/summary.json)

The frozen circuit does not reproduce the complex HFSS S21. The deepest notch is 42 MHz below HFSS and the linear complex RMS difference is 0.354 over 6.70–7.08 GHz. The direct comparison uses matching signal-geometry snapshots and ground dimensions with the fixed port polarity. [Amplitude and phase overlay](results/resonator/20261118_170023_q3d_circuit_comparison_a4edb91d/s21_q3d_vs_hfss.png) · [Run index and raw exports](results/README.md)

[Reproduction in qmetal](simulation/resonator/README.md)
