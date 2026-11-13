# Design3 現版電性確認 · 2026-11-13

下 pad 耦合，paddle gap 20 µm，TC readout pad 150 × 20 µm，矩形 SQUID 開孔 20 × 14 µm。矽基板模型為 525 µm。

## Q3D 電容與 readout 耦合

| 裝置 | EC/h (GHz) | CΣ (fF) | g/2π 每單位 n (MHz) | g01/2π (MHz) |
|---|---:|---:|---:|---:|
| FX1 | 1.0623 | 18.234 | 87.615 | 8.475 |
| TC | 0.2149 | 90.148 | 114.382 | 137.921 |
| FX2 | 1.0643 | 18.201 | 87.356 | 8.466 |

TC 依 Design3 原本的 EJ/EC = 50 慣例為 **102.306 MHz**，目標為 100 MHz。
g01 採 FX 半磁通、EJ/h = 4 GHz、EL/h = 0.6 GHz，TC f01 = 5 GHz。readout 模型採半波 CPW 與 50.9 Ω。頻率取下表的 readout 模式估值。

| Qubit pair | J/2π (MHz) | 0→1 耦合幅度 (MHz) |
|---|---:|---:|
| FX1-TC | +127.3563 | 14.8551 |
| FX2-TC | -125.7110 | 14.6899 |
| FX1-FX2 | -9.3164 | 0.0873 |

### Charge-line 耦合

| 接收端 | FX1 line β | TC line β | FX2 line β |
|---|---:|---:|---:|
| FX1 | 0.007641 | 2.546e-05 | 8.3925e-06 |
| TC | 0.00065093 | 0.0003309 | 0.00066114 |
| FX2 | 3.5355e-06 | -0.00017334 | -0.0078933 |

三條 charge line 對 TC 的直接電容依 FX1、TC、FX2 排列為 0.05927 fF / 0.03013 fF / 0.06020 fF。

## HFSS readout 頻率

| 裝置 | 目標 (GHz) | readout 估值 (GHz) | 偏差 (MHz) | 最後一 pass 漂移 (MHz) | HFSS 收斂 | 全部候選 (GHz) |
|---|---:|---:|---:|---:|---|---|
| FX1 | 6.300 | 6.279244 | -20.76 | +7.55 | 未達門檻 | 6.279244, 6.498317 |
| TC | 6.000 | 6.057906 | +57.91 | -6.11 | 未達門檻 | 6.057906, 8.712198 |
| FX2 | 6.500 | 6.457098 | -42.90 | +9.01 | 未達門檻 | 6.457098, 6.573284 |

各 resonator 個別求解，包含其 tee 和 qubit，TC 兩顆 JJ 各用 50 nH。頻率為線性模型、有限 pass 的估值。FX1、FX2 採較低的候選模式，兩者在 CPW 槽附近的電場取樣均較強。TC 取接近半波設計的低頻候選，兩個 TC 模式的取樣皆顯示 CPW 場，模式辨認仍有限。

## Resonator–feedline 耦合

| 裝置 | 耦合電容 (fF) | κ/2π 估值 (MHz) |
|---|---:|---:|
| FX1 | 15.6457 | 1.9731 |
| TC | 20.5280 | 3.0500 |
| FX2 | 14.8599 | 1.9355 |

κ 由被動 feedline Q3D 模型和半波線公式估算。12 pass 時矩陣 ΔC = 1.5293%，未達 0.2% 設定。

Q3D 矩陣在 6 pass 停止，最後 ΔC = 4.3501%，未達 0.5%。第 5→6 pass，EC/h 最大差 14.05 MHz，J 最大差 3.28 MHz。SQUID 接地片與地平面在 AEDT 先合併成同一導體。

## 原始資料

- [完整數據及 3 × 3 矩陣](simulation/results/runs/d3_final_parameters.json)
- [Maxwell 電容矩陣](simulation/results/runs/d3_final_matrix_cmat_fF.csv)
- [Q3D 收斂](simulation/results/runs/d3_final_matrix_convergence.csv)
- [HFSS 專案](simulation/aedt/Design3.aedt)
- [先前版圖與互感比較](design_change_comparison.md)
- FX1 [頻率及模式](simulation/results/runs/d3_final_eigen_FX1.csv) · [pass 曲線](simulation/results/runs/d3_final_eigen_FX1_convergence.csv) · [場取樣](simulation/results/runs/d3_final_eigen_FX1_field_samples.csv)
- TC [頻率及模式](simulation/results/runs/d3_final_eigen_TC.csv) · [pass 曲線](simulation/results/runs/d3_final_eigen_TC_convergence.csv) · [場取樣](simulation/results/runs/d3_final_eigen_TC_field_samples.csv)
- FX2 [頻率及模式](simulation/results/runs/d3_final_eigen_FX2.csv) · [pass 曲線](simulation/results/runs/d3_final_eigen_FX2_convergence.csv) · [場取樣](simulation/results/runs/d3_final_eigen_FX2_field_samples.csv)
