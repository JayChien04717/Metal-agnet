# r9 TC mutual inductance

The r9 local TC fixture gives mutual-inductance magnitudes of 1.189439 pH at the left JJ sense gap and 1.161899 pH at the right gap at 10 MHz. Both exceed the requested 1 pH threshold. The mean magnitude is 1.175669 pH.

The model uses the imported TC electrodes, the 3.5 µm displaced flux termination and an explicit PEC closure across the other JJ. The drive port is at the local straight feed near (-55, -1800) µm. The shorted flux conductor and ground are united in HFSS. The simulation uses the existing 10–100 MHz discrete sweep and 15-pass maximum, with 12 minimum adaptive passes. The local left fixture completed 12 passes with final ΔS = 0.000064204. Pass-count stability beyond 12 passes has not been measured.

[Recorded run](../results/flux/20261118_181301_r9_mutual_027c2275/run.json) · [Summary](../results/flux/20261118_181301_r9_mutual_027c2275/summary.json) · Left raw data (`results/flux/20261118_181301_r9_mutual_027c2275/raw/D2_M_tc_left.csv`) · Right raw data (`results/flux/20261118_181301_r9_mutual_027c2275/raw/D2_M_tc_right.csv`) · [Reproduction](../simulation/flux/recheck_ground.py)

Reproduce with `python -B Design2_FTF_Andersen/simulation/flux/recheck_ground.py --local` in qmetal. The local model excludes the long route and launch. It checks the device-end geometry.

Two diagnostic runs used the legacy full-route fixture. The first had a launch port sheet ending before ground. Extending that port restored inductive self-impedance, but the fixture still mixed the original conductor with a rebuilt launch ground. Those values are unsuitable for adopting the current wafer's mutual inductance. Their immutable runs remain at `results/flux/20261118_175435_r9_mutual_11a63db0` and `results/flux/20261118_180833_r9_mutual_deb31611`. The local run above is the selected device-end check.
