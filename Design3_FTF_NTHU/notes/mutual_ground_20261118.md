# r9 TC mutual inductance

The r9 geometry retains the r8 TC flux line with 10.5 µm separation. Local HFSS fixtures use one JJ gap as a sense port and close the other gap with PEC. Both cases retain the existing 10–100 MHz discrete sweep and 0.0002 adaptive ΔS target.

| Sense gap | Passes completed | Mutual-inductance magnitude at 10 MHz (pH) | Final ΔS |
|---|---:|---:|---:|
| Left | 12 | 0.988036 | 0.000016781 |
| Right | 12 | 1.068375 | 0.000019627 |
| Left, finer mesh | 15 | 1.002645 | 0.0000099417 |

The finer left result exceeds 1 pH by about 0.26%. The 12-to-15-pass change is about 1.48%, so the available mesh evidence does not establish a margin above 1 pH. The left and right fixtures also give different values because the port and PEC closure exchange sides. The original report's mean-magnitude convention should be stated explicitly when quoting one TC value.

[12-pass run](../results/flux/20261118_174738_r9_mutual_71b19ea5/run.json) · [15-pass run](../results/flux/20261118_175915_r9_mutual_326b3e4b/run.json) · 15-pass raw data (`results/flux/20261118_175915_r9_mutual_326b3e4b/raw/m_tc_left_Z.csv`) · [Reproduction](../simulation/flux/recheck_ground.py)

The simulations are local TC fixtures. Full-chip routing and fabrication tolerance are outside this check.
