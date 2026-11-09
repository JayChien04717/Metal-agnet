# pattern design（共用 QComponent）

這裡只放元件（QComponent）和 port 位置，供 Design1–6 共用。晶片 layout 留在各 Design 自己的
`pattern design/`（`ftf_layout.py`、`purcell_layout.py`、`design5_layout.py`、`chip6_layout.py`…），
那些 layout 會把本資料夾加進 sys.path。

| 類別 | 檔案 | 類別名 |
|---|---|---|
| port | `pcbport.py` | `pcb_ports`（PCB port 位置） |
| fluxonium | `my_fluxonium.py` | `MyFluxonium`（預設 flux line 80 / 11 µm，Design3；Design1 在 `FX_OPTIONS` 寫 50 / 15 µm） |
| | `fluxonium.py` | `FluxoniumArrayJJ`（量互感用） |
| | `fluxonium_forprocess.py` | `FluxoniumForProcess`（製程版） |
| | `fluxonium_withJJ.py` | `FluxoniumPocket`（Design2） |
| transmon / coupler | `my_tmon.py` | `MyTmonCoupler` |
| | `tmon_c.py` | `Tmon_c`（Design2） |
| 電容 | `hammer_cap.py` | `HammerCap`（C_j） |
| | `enclosing_cap.py` | `EnclosingCap`（hugger，C_in / C_p） |
| | `interdigital_cap.py` | `InterdigitalCap` |
| 其他 | `tap_junction.py`, `res_end.py`, `test_junction.py` | T 接頭、開路 / 短路端、測試 junction |
