# Design10_CPW_BPF_5p5_6p5GHz Results

| Run ID | Type | Status | Notes |
|---|---|---|---|
| [`20261121_162252_hfss_driven_1port_3stage_9d5ee50f`](resonator/20261121_162252_hfss_driven_1port_3stage_9d5ee50f) | HFSS Driven (`2.0–7.5 GHz`, Fast sweep) | Raw complete (`summaries/` missing at post-process) | 3 passes (user override), `10 µm` mesh, 40 PEC ground wirebonds, `451,272` tets, Pass 3 `Max Mag. Delta S = 0.033366`. Contains immutable `raw/S_params_raw50.csv` and `raw/hfss_convergence.csv`. |
| [`20261121_171558_reanalysis_hfss_1port_s21_s11_c4dd74ce`](resonator/20261121_171558_reanalysis_hfss_1port_s21_s11_c4dd74ce) | Re-analysis of parent HFSS run | **Adopted (`completed`)** | Clean pure-HFSS exports: `summaries/hfss_1port_2_7p5GHz.csv`, `summaries/hfss_1port_s11.s1p`, `figures/hfss_1port_s21_s11.png`. `-3 dB`: `5.358–6.560 GHz`, `5.5–6.5 GHz` min `-1.768 dB`, `2–4 GHz` max `-43.56 dB`. |
