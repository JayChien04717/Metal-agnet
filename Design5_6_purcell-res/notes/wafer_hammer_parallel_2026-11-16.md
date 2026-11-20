# Design5_6_purcell-res straight parallel CPW hammer

The r5 wafer candidate uses two straight 16/9 µm CPWs at a centre spacing
of 190 µm through each hammer. The inner etched edges are 156 µm apart.
For each of the six channels the filter stays on its original axis. The
readout resonator and its open end and meander translate 210 µm toward
the filter. The existing filter and resonator centreline lengths are
retained. This removes the local S bends and their extra path length.
The intended coupling pair uses the user-specified 190 µm separation.
The layout check measured a minimum 878.9 µm metal separation between
different channels.

The local Q3D fixture has 439.6 µm straight CPWs on the same substrate
and four open ends. Its 400 µm reference geometry matches the saved
reference exactly and is reused from the [r4 Q3D run](../results/resonator/20261116_160806_hammer_190_q3d_070bbb59/run.json).
The straight comparison uses the same CPW lengths and extraction settings.
The selected head is 25 × 52.2 µm, with a 40 µm face gap and 50 µm from
each CPW centre to the outside of its head. Its mutual capacitance is
1.50725 fF against 1.50888 fF for the original straight 400 µm reference,
a difference of −0.11%. The selected case reports solver convergence.
The last two mutual estimates differ by 0.76% for the selected case and
0.83% for the reference. The reference did not report convergence.
See the [Q3D summary](../results/resonator/20261116_161734_hammer_parallel_190_q3d_dd191dce/summary.json)
and [run manifest](../results/resonator/20261116_161734_hammer_parallel_190_q3d_dd191dce/run.json).

The final check reads back all six CPW pairs from the delivered mask and models
their actual 1,317 µm shared straight span. The mutual capacitance is 2.06857 fF
for r5 and 1.68841 fF for the straight 400 µm reference, a +22.52% difference.
Both simulations saved eight adaptive passes. The last two mutual estimates
differ by 0.37% for r5 and 0.68% for the reference. Neither reports solver
convergence. The current [r7 fab mask](../gds/Design5_2026-11-16_wafer_r7_fab.gds) uses the same
readout geometry and the 20 µm-pitch ground perforations. See the
[final Q3D comparison](../results/resonator/20261116_163137_hammer_final_overlap_q3d_f73e11bd/summary.json).

Both comparisons omit the rest of the filter and resonator and the ground
bridges. Full-chip coupling and resonant frequencies still require electrical
revalidation. HFSS was not rerun.

The fab GDS uses layer 1/0 positive metal, µm units and a 10 × 10 mm
outline. It has centred launches and the standard ground perforation.
It omits the cross markers, die ID and corner test structures.

Run from the repository root in the `qmetal` environment:

```powershell
python -B Design5_6_purcell-res/simulation/resonator/q3d_hammer_check.py
python -B Design5_6_purcell-res/simulation/resonator/q3d_hammer_final.py
python -B code/general/export_wafer_test_junctions.py --design Design5_6_purcell-res --output Design5_r7_fab_check.gds
```

The wafer export reads the selected geometry from [design.json](../design.json).
