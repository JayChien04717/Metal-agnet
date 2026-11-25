# E-beam junction masks

These are the shared QComponents for the Dolan-bridge junction layers. The rules, the EJ → resistance → area → width
model, the area-vs-resistance table and the workflow are in [skill/ebeam_junction.md](../../skill/ebeam_junction.md).

| Component | Draws | Options |
|---|---|---|
| `FluxoniumJunctionDolan` ([py](fluxonium_junction_dolan.py)) | Fluxonium small JJ, JJ array, funnel contacts | [fluxonium_junction_dolan.md](fluxonium_junction_dolan.md) |
| `TransmonJunctionDolan` ([py](transmon_junction_dolan.py)) | Transmon / coupler Dolan SQUID or single JJ | [transmon_junction_dolan.md](transmon_junction_dolan.md) |
| `TransmonJunctionManhattan` ([py](transmon_junction_manhattan.py)) | Transmon / coupler Manhattan SQUID or single JJ (45° pads + orthogonal fingers + undercut) | [transmon_junction_manhattan.md](transmon_junction_manhattan.md) |

Sizing comes from [`code/qubit/jj_resistance_model.py`](../../code/qubit/jj_resistance_model.py), and the checks from
[`code/ebeam/check_ebeam_junctions.py`](../../code/ebeam/check_ebeam_junctions.py).
