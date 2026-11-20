# E-beam junction masks

These are the shared QComponents for the Dolan-bridge junction layers. The rules, the EJ → resistance → area → width
model, the area-vs-resistance table and the workflow are in [skill/ebeam_junction.md](../../skill/ebeam_junction.md).

| Component | Draws | Options |
|---|---|---|
| `FluxoniumJunction` ([py](fluxonium_junction.py)) | Fluxonium small JJ, JJ array, funnel contacts | [fluxonium_junction.md](fluxonium_junction.md) |
| `TransmonJunction` ([py](transmon_junction.py)) | Transmon / coupler SQUID or single JJ | [transmon_junction.md](transmon_junction.md) |

Sizing comes from [`code/qubit/jj_resistance_model.py`](../../code/qubit/jj_resistance_model.py), and the checks from
[`code/general/check_ebeam_junctions.py`](../../code/general/check_ebeam_junctions.py). The old layout3 drawings are in
[`../archive/`](../archive/README.md).
