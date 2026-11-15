# Design3_2 模擬程式

按研究目的分類，不按 solver 分類。general 是跨用途檢查；HFSS／Q3D／scqubits 名稱保留在入口檔名中。

| 入口 | 分類 |
|---|---|
| [prepare_fx1_port_review.py](flux/prepare_fx1_port_review.py) | flux |
| [run_remaining.py](general/run_remaining.py) | general |
| [run_targeted.py](general/run_targeted.py) | general |
| [run_verification.py](general/run_verification.py) | general |
| [update_readme.py](general/update_readme.py) | general |

舊入口的輸入／輸出介面保留；部分腳本仍會寫固定檔名，執行前需指定新的輸出位置。`code/general/run_record.py` 可建立新 run 並封存程式／參數，但不會暗中重導舊腳本輸出。來源分類只調整路徑，不表示重跑或驗證結果。
