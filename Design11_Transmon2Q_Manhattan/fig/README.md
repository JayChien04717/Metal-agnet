Store full-chip/local figures for the same revision, following skill/layout_figures.md. Do not add placeholder images before generation.

Current files:
- Layout figures (`r17`, not simulated, rendered by `layout/build.py` from `gds/Design11_2026-11-26_r17_sim.gds`): `chip.png` (whole die, 8 ports), `core.png` (window x ±0.75, y ±0.45 mm), `squid_q1.png` (x -0.45 … -0.35, y -0.24 … -0.14 mm). The `r9` and `r15` versions are in `previous/` (ignored).
- HFSS readout spectrum (`r15` geometry, accepted 2026-11-26):
  - `readout_s21_r15.png`: full-band `|S21|` / `|S11|`. Plotted from `results/resonator/20261126_074508_purcell_s21_driven_r15_kp30_3p_fine4_afc42719/raw/S_wide.csv`.
  - `readout_s21_fit_r15.png`: complex fit. Copied from `results/resonator/20261126_081341_fit_s21_r15_kp30_3p_fine4_9b88474c/figures/s21_fit.png`.
- CZ simulation on the `r9` core (unchanged in `r15`), copied from `results/qubit/20261125_220541_cz_verify_r9_optB_289426f4/figures/`:
  - `cz_vs_time_r9core.png`
  - `zz_vs_flux_r9core.png`
