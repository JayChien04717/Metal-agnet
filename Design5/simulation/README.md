# Design5 模擬程式

按研究目的分類，不按 solver 分類。general 是跨用途檢查；HFSS／Q3D／scqubits 名稱保留在入口檔名中。

| 入口 | 分類 |
|---|---|
| [analyze5.py](resonator/analyze5.py) | resonator |
| [chip6_d4_full.py](resonator/chip6_d4_full.py) | resonator |
| [chip6_kappa_env.py](resonator/chip6_kappa_env.py) | resonator |
| [chip6_theory.py](resonator/chip6_theory.py) | resonator |
| [eig_chip6.py](resonator/eig_chip6.py) | resonator |
| [hfss_chan5.py](resonator/hfss_chan5.py) | resonator |
| [hfss_chip6.py](resonator/hfss_chip6.py) | resonator |
| [modes5.py](resonator/modes5.py) | resonator |
| [q3d_full5.py](resonator/q3d_full5.py) | resonator |
| [runstore.py](resonator/runstore.py) | resonator |
| [solve_chip6.py](resonator/solve_chip6.py) | resonator |
| [theory5.py](resonator/theory5.py) | resonator |

舊入口的輸入／輸出介面保留；部分腳本仍會寫固定檔名，執行前需指定新的輸出位置。`code/general/run_record.py` 可建立新 run 並封存程式／參數，但不會暗中重導舊腳本輸出。來源分類只調整路徑，不表示重跑或驗證結果。
