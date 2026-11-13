# pattern design（共用 QComponent）

這裡只放元件（QComponent）和 port 位置，供 Design1–6 共用。晶片 layout 留在各 Design 自己的
`pattern design/`（`ftf_layout.py`、`purcell_layout.py`、`design5_layout.py`、`chip6_layout.py`…），
那些 layout 會把本資料夾加進 sys.path。

| 類別 | 檔案 | 類別名 |
|---|---|---|
| port | `pcbport.py` | `pcb_ports`（PCB port 位置） |
| fluxonium | `my_fluxonium.py` | `MyFluxonium`（預設 flux line 80 / 11 µm，Design3；Design1 在 `FX_OPTIONS` 寫 50 / 15 µm） |
| | `fluxonium.py` | `FluxoniumArrayJJ`（所有 fluxonium 模擬用，loop 同 `MyFluxonium._loop`：JJ = lead 線中間一小段（qubit 中心），U 的兩側是 lumped L（各 100 nH），U 末端是金屬；Q3D 用 `flux_as_ground='True'`） |
| | `fluxonium_withJJ.py` | `FluxoniumPocket`（Design2；模擬用 `u_array='True'`：同樣的 loop） |
| transmon / coupler | `my_tmon.py` | `MyTmonCoupler` |
| | `tmon_c.py` | `Tmon_c`（Design2） |
| 電容 | `hammer_cap.py` | `HammerCap`（C_j） |
| | `enclosing_cap.py` | `EnclosingCap`（hugger，C_in / C_p） |
| | `interdigital_cap.py` | `InterdigitalCap` |
| 其他 | `tap_junction.py`, `res_end.py`, `test_junction.py` | T 接頭、開路 / 短路端、測試 junction |
| GDS | `gds_tools.py` | `export_gds`（junction-free）、`plot_gds`、`junction_window`（Design1–3 共用，圖格式一致；見 `skill/gds_export.md`、`skill/layout_figures.md`） |
| 標記 | `xmark.py` | `XMark`（蝕刻正方形窗口 `window_size` 400 µm，中間是金屬十字 250 × 10 µm；`label` 蝕刻在窗口右邊、十字線下方，間距 `label_gap` 30 µm，字距緊湊）、`add_corner_marks()`（晶片四角各一個，中心距邊 0.35 mm，die ID 放在左上角那個） |
