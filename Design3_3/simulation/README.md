# Design3_3 模擬程式

按研究目的分類，不按 solver 分類。general 是跨用途檢查；HFSS／Q3D／scqubits 名稱保留在入口檔名中。

| 入口 | 分類 |
|---|---|
| [recheck_readout_flux.py](flux/recheck_readout_flux.py) | flux |
| [report_readout_flux.py](flux/report_readout_flux.py) | flux |
| [verify.py](general/verify.py) | general |
| [compare.py](qubit/compare.py) | qubit |
| [optimize_pulse.py](qubit/optimize_pulse.py) | qubit |
| [q3d.py](qubit/q3d.py) | qubit |
| [recheck_pulse.py](qubit/recheck_pulse.py) | qubit |
| [run_q3d.py](qubit/run_q3d.py) | qubit |
| [sequence.py](qubit/sequence.py) | qubit |
| [sweep.py](qubit/sweep.py) | qubit |
| [validate_gate.py](qubit/validate_gate.py) | qubit |
| [probe_readout.py](resonator/probe_readout.py) | resonator |

舊入口的輸入／輸出介面保留；部分腳本仍會寫固定檔名，執行前需指定新的輸出位置。`code/general/run_record.py` 可建立新 run 並封存程式／參數，但不會暗中重導舊腳本輸出。來源分類只調整路徑，不表示重跑或驗證結果。
