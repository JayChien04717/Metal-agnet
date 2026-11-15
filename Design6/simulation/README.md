# Design6 模擬程式

按研究目的分類，不按 solver 分類。general 是跨用途檢查；HFSS／Q3D／scqubits 名稱保留在入口檔名中。

| 入口 | 分類 |
|---|---|
| [compare_flux_variants.py](flux/compare_flux_variants.py) | flux |
| [flux_metrics.py](flux/flux_metrics.py) | flux |
| [flux_pulse.py](flux/flux_pulse.py) | flux |
| [flux_resweep.py](flux/flux_resweep.py) | flux |
| [flux_square_spectro.py](flux/flux_square_spectro.py) | flux |
| [flux_tail_scan.py](flux/flux_tail_scan.py) | flux |
| [flux_xtalk_matrix.py](flux/flux_xtalk_matrix.py) | flux |
| [hfss_flux.py](flux/hfss_flux.py) | flux |
| [cz_scqubits.py](qubit/cz_scqubits.py) | qubit |
| [hfss_core.py](qubit/hfss_core.py) | qubit |
| [hfss_pad.py](qubit/hfss_pad.py) | qubit |
| [plot_ej_layout.py](qubit/plot_ej_layout.py) | qubit |
| [q3d_core.py](qubit/q3d_core.py) | qubit |
| [epr_readout.py](resonator/epr_readout.py) | resonator |
| [hfss_readout.py](resonator/hfss_readout.py) | resonator |

舊入口的輸入／輸出介面保留；部分腳本仍會寫固定檔名，執行前需指定新的輸出位置。`code/general/run_record.py` 可建立新 run 並封存程式／參數，但不會暗中重導舊腳本輸出。來源分類只調整路徑，不表示重跑或驗證結果。
