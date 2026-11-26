# Design12 GDS (local only, not on GitHub)

| File | Role | Made by |
|---|---|---|
| `Design12_Transmon2Q_TransmissionPurcell_2026-11-26_r4_sim.gds` | simulation reference of the adopted `r4` (layer 1, lumped junction sheets on 50/0) | `layout/export_fab.py` |
| `Design12_Transmon2Q_TransmissionPurcell_2026-11-26_r4_fab.gds` | NCKU single-die fab mask, die ID `prt2q00`: `jj_metal=False` (junction sites empty), R = 5 um rounding, 5 x 5 um / 40 um cheese, no markers, no test junctions | `layout/export_fab.py` |
| `previous/` | superseded `r1`-`r3` sim GDS | `layout/build.py` |

Checks and SHA256 of the current pair: the latest `results/general/*_gds_export_r4_ncku_*` run.
