# r11 transmon mutual inductance

The confirmed r11 flux-line geometry uses 10.5 µm separation and 6 µm CPW slots. The local three-transmon HFSS fixture gives the following magnitudes at 10 MHz.

| Device | Mutual-inductance magnitude (pH) | Current per flux quantum (mA) |
|---|---:|---:|
| Q1 | 1.17587 | 1.75856 |
| TC | 1.24428 | 1.66187 |
| Q2 | 1.17121 | 1.76556 |

All three magnitudes exceed the requested 1 pH threshold. The signed values depend on port direction. Six adaptive passes with a 4 µm seed mesh ended at ΔS = 0.000023155 against the 0.0002 setup target. A second pass-count comparison has not been performed.

The fixture includes Q1, TC and Q2 and their local flux lines, with ports at the local pins. It omits the long routed bias lines. The discrete sweep spans 10–1000 MHz. Its full-band variation includes higher-frequency effects and is recorded in the raw matrix summary.

Raw solve (`results/flux/20261118_171936_r11_mutual_9d9821fe/raw/M_matrix.json`) · Convergence (`results/flux/20261118_171936_r11_mutual_9d9821fe/raw/convergence.csv`) · [Recovered summary](../results/flux/20261118_174849_r11_mutual_summary_b46b1b70/summary.json) · [Reproduction](../simulation/flux/recheck_ground.py)

The solver finished and saved its AEDT snapshot. The initial summary writer failed on a NumPy boolean after storing the raw data. A separate summary run records the recovered result and its parent solve.
