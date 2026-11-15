# Design2 模擬程式

按研究目的分類，不按 solver 分類。general 是跨用途檢查；HFSS／Q3D／scqubits 名稱保留在入口檔名中。

| 入口 | 分類 |
|---|---|
| [hfss_mutual_fx.py](flux/hfss_mutual_fx.py) | flux |
| [run_mutual.py](flux/run_mutual.py) | flux |
| [analyse.py](general/analyse.py) | general |
| [q3d_ftf.py](qubit/q3d_ftf.py) | qubit |
| [run_q3d.py](qubit/run_q3d.py) | qubit |
| [transition_estimates.py](qubit/transition_estimates.py) | qubit |
| [hfss_eigen_res.py](resonator/hfss_eigen_res.py) | resonator |
| [linewidth.py](resonator/linewidth.py) | resonator |
| [modal_linewidth.py](resonator/modal_linewidth.py) | resonator |
| [plot_readout_modes.py](resonator/plot_readout_modes.py) | resonator |
| [q3d_kappa.py](resonator/q3d_kappa.py) | resonator |
| [readout_local_eigen.py](resonator/readout_local_eigen.py) | resonator |
| [run_linear_readout.py](resonator/run_linear_readout.py) | resonator |
| [run_readout.py](resonator/run_readout.py) | resonator |
| [run_res_modes.py](resonator/run_res_modes.py) | resonator |

舊入口的輸入／輸出介面保留；部分腳本仍會寫固定檔名，執行前需指定新的輸出位置。`code/general/run_record.py` 可建立新 run 並封存程式／參數，但不會暗中重導舊腳本輸出。來源分類只調整路徑，不表示重跑或驗證結果。
