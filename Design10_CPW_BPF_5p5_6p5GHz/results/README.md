# Design10_CPW_BPF_5p5_6p5GHz Results

| Run ID | Type | Status | Notes |
|---|---|---|---|
| [`20261121_162252_hfss_driven_1port_3stage_9d5ee50f`](resonator/20261121_162252_hfss_driven_1port_3stage_9d5ee50f) | HFSS Driven (`2.0–7.5 GHz`, Fast sweep) | Raw complete (`summaries/` missing at post-process) | 3 passes (user override), `10 µm` mesh, 40 PEC ground wirebonds, `451,272` tets, Pass 3 `Max Mag. Delta S = 0.033366`. Contains immutable `raw/S_params_raw50.csv` and `raw/hfss_convergence.csv`. |
| [`20261121_171558_reanalysis_hfss_1port_s21_s11_c4dd74ce`](resonator/20261121_171558_reanalysis_hfss_1port_s21_s11_c4dd74ce) | Re-analysis (50 kΩ probe termination) | Superseded by ideal-open plot | Prior 50 kΩ probe termination (`|S11|` had `-0.15 dB` probe loading). |
| [`20261121_210808_reanalysis_hfss_1port_s11_vopen_c8ba26ca`](resonator/20261121_210808_reanalysis_hfss_1port_s11_vopen_c8ba26ca) | Re-analysis (ideal open $\Gamma_2 = +1$) | **Adopted (`completed`)** | Circulator+VNA `S21 = S11` (`0.000 dB`, phase, $\tau_g$) + internal open-end `|V_open/V_inc|` (`-3 dB`: `5.360–6.558 GHz`, `5.5–6.5 GHz` min `-1.798 dB`, `2–4 GHz` max `-43.65 dB`). |
