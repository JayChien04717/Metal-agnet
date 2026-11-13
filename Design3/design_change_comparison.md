# Design3 版圖與互感比較 · 2026-11-12

[目前參數](design_parameters.txt) · [原參數檔快照](simulation/results/previous/design3_parameters_before_rect_squid.txt) · [HFSS 專案](simulation/aedt/Design3.aedt) · [重現腳本](simulation/driven/run_mutual_current.py)


## Junction-free GDS

- [一般版](gds/FTF00_metal_only.gds) · [有地孔的 fab 版](gds/FTF00_metal_only_fab.gds)
- [匯出前的一般版](gds/previous/ftf_layout.gds) · [匯出前的 fab 版](gds/previous/ftf_layout_fab.gds)
- 目前 FTF00 版保留 TC SQUID 底層落腳電極，兩個 5 µm JJ 間隙留空。JJ 圖形與 fluxonium array 由製程方補畫。SQUID 朝 stem 移 2 µm，接地端重疊為 1 µm。
## 版圖變更

| 參數 | 舊版 | 新版 | 差異 |
|---|---|---|---|
| TC paddle 對應 FX 電極 | 上 pad | 下 pad | bar 中心 y 從 −0.715 移到 −0.885 mm |
| TC–FX paddle 金屬間隙 | 20 µm | 20 µm | 0 µm |
| FX readout 位置 | 上 pad | 上 pad | 0 µm |
| TC readout pad | 273 × 20 µm | 150 × 20 µm | 寬度 −123 µm（−45.1%），gap 5 µm |
| TC SQUID 開孔 | 16 × 24 µm 外接六角孔，288 µm² | 20 × 14 µm 矩形，280 µm² | −8 µm²（−2.8%） |
| SQUID 電極和底層金屬 | 邊對邊接觸 | island、ground 端各 3 µm 名義重疊 | +3 µm |
| TC SQUID 兩顆 JJ 的中心間距 | 28 µm | 30 µm | +2 µm |
| FX flux line | U 腿間距 80 µm，tip gap 11 µm | U 腿間距 80 µm，tip gap 11 µm | 0 µm |

TC 100 MHz 是 readout 耦合目標。下 pad 與縮小 readout pad 後的電容耦合仍待 Q3D 求解。舊版電容數據保留在 [上 pad 矩陣](simulation/results/runs/d3_pw152_ro21_cmat_fF.csv) 和 [原摘要](simulation/results/runs/d3_pw152_ro21_summary.json)。

## 自身 bias line → loop 互感

| 裝置 | 舊 | 新 | 差值 |
|---|---:|---:|---:|
| TC，兩版均用左 JJ 開口 | 0.8026 pH | 0.7322 pH | -0.0704 pH（-8.8%） |
| FX1 | 1.0369 pH | 1.0307 pH | -0.0062 pH（-0.6%） |
| FX2 | 無獨立歷史結果 | 1.0313 pH | — |

TC 矩形左、右兩個取樣的幅度平均 **0.7175 pH**，兩次全差 **4.1%**。以平均值換算，半磁通需 **1.4410 mA**。六角孔基準是在目前其他 TC 幾何上還原舊開孔與接觸設定後重跑。FX1 舊值來自歷史設計，版圖改動與元件建模差異都可能影響差值。

## 原始曲線與收斂

| Case | 10 MHz 帶符號 M（pH） | 半磁通（mA） | 10–100 MHz M 變幅 | Passes | 最後 ΔS | 資料 |
|---|---:|---:|---:|---:|---:|---|
| TC old | -0.8026 | 1.2881 | 0.0001722% | 12 | 3.15e-05 | [Z](simulation/results/runs/d3_tc_hex16x24_m_left_Z.csv) · [summary](simulation/results/runs/d3_tc_hex16x24_m_left_summary.json) · [mesh](simulation/results/runs/d3_tc_hex16x24_m_left_convergence.csv) |
| TC left | -0.7322 | 1.4120 | 0.000139% | 12 | 4.1e-05 | [Z](simulation/results/runs/d3_tc_rect20x14_m_left_Z.csv) · [summary](simulation/results/runs/d3_tc_rect20x14_m_left_summary.json) · [mesh](simulation/results/runs/d3_tc_rect20x14_m_left_convergence.csv) |
| TC right | +0.7028 | 1.4712 | 0.0004564% | 12 | 3.85e-05 | [Z](simulation/results/runs/d3_tc_rect20x14_m_right_Z.csv) · [summary](simulation/results/runs/d3_tc_rect20x14_m_right_summary.json) · [mesh](simulation/results/runs/d3_tc_rect20x14_m_right_convergence.csv) |
| FX1 old | +1.0369 | 0.9971 | 1.124% | 12 | 3.43e-07 | [Z](simulation/results/runs/mutual_d3_fl11_h80_Z.csv) · [summary](simulation/results/runs/mutual_d3_fl11_h80_summary.json) · [mesh](simulation/results/runs/mutual_d3_fl11_h80_convergence.csv) |
| FX1 | +1.0307 | 1.0031 | 1.4% | 12 | 2.59e-07 | [Z](simulation/results/runs/d3_current_m_fx1_Z.csv) · [summary](simulation/results/runs/d3_current_m_fx1_summary.json) · [mesh](simulation/results/runs/d3_current_m_fx1_convergence.csv) |
| FX2 | -1.0313 | 1.0025 | 1.361% | 12 | 1.33e-06 | [Z](simulation/results/runs/d3_current_m_fx2_Z.csv) · [summary](simulation/results/runs/d3_current_m_fx2_summary.json) · [mesh](simulation/results/runs/d3_current_m_fx2_convergence.csv) |

HFSS Driven Modal 在 100 MHz 自適應，10–100 MHz 做 10 點離散掃頻。以 (M=\mathrm{Im}[(Z_{21}+Z_{12})/2]/(2\pi f)) 取 10 MHz 值。TC 將一側 JJ 設為開路量測端，另一側以 PEC 封閉。FX 保留兩個 100 nH array 臂與金屬封口，主 JJ 為開路量測端。TC 左、右開口的符號相反，反映端口的取向。模型採單元附近的 ground return。結果適用於自身 flux line 的互感。
