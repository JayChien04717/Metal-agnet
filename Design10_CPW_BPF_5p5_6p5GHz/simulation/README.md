# Design10_CPW_BPF_5p5_6p5GHz Simulation

- **HFSS Driven 1-port simulation (`2.0–7.5 GHz`, Fast sweep):** [`resonator/hfss_s21.py`](resonator/hfss_s21.py)
  - User-requested override: `max_passes = 3` (repo default is 4) after explicit `10 µm` mesh on all 10 signal objects + 40 PEC ground wirebonds.
  - Physical Port 1: `('feed_in', 'start', 50)`. Non-perturbing voltage probe at `R3_open`: `('R3', 'end', 50000)`.
- **Pure HFSS post-process / re-analysis exporter (CSV, `.s1p`, HFSS-only plot):** [`resonator/export_hfss_1port.py`](resonator/export_hfss_1port.py)
