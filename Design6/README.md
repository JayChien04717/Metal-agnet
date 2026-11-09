# Design6：驗證電容模型能不能對上 HFSS S21 的耦合

（2026-11-09 起。舊的 Design6（Cai 2024 型 channel、c6_* runs）已刪除，需要時從 git 歷史取回。）

## 要回答的問題

1. **Purcell ↔ feed line**：Q3D 算出的 C_p，放進電路模型後預測的 ω_p、κ_p，和 HFSS S21 一致嗎？
2. **resonator ↔ Purcell**：Q3D 算出的 C_j，預測的 J，和 HFSS S21 一致嗎？

反過來說：把電路模型擬合到 HFSS S21，得到的 C_p、C_j 和 Q3D 的值差多少？

## 為什麼這樣設計

- **沒有 C_in**：feed line 從左邊 pad 直接連到右邊 pad，兩端匹配。S21 = 2 / (2 + Z0·Y_branch)（只差一個相位），
  κ_p 不會被駐波改變（有 C_in 時 κ_p 隨位置 0.9–93 MHz，見 Design5 README 9b）。
- **一次只驗證一個耦合**：
  - `filter`：只有 C_p + Purcell filter（resonator 的線拿掉，C_j hammer 保留、另一片浮接，所以 filter 的負載和完整版一樣）。
    S21 只有**一個** notch → ω_p、κ_p 直接看得到，只和 C_p（和 filter 長度）有關。
  - `full`：完整 channel。C_p 已在上一步驗證，剩下 C_j → J。
- 單一 channel（6.7 GHz，Design5 chip6 ch3 的長度），hairpin 走線。

## 檔案

```
Design6/
├─ README.md
├─ pattern design/verify_layout.py   layout：VARIANT filter / full（env CP_PLATE、CJ_GAP 可做電容掃描）
├─ simulation/verify_model.py        電路模型（沒有 C_in）與由電容直接得到的預測
├─ gds/  verify_filter.gds, verify_full.gds
└─ fig/  verify_filter(_zoom).png, verify_full(_zoom).png
```

## 步驟

| 步驟 | 內容 | 比較的東西 |
|---|---|---|
| 1. Q3D | `verify_full` 整條 channel 一起跑 | 這個幾何（沒有 C_in）的 C_p、C_j 和對地 C |
| 2. HFSS `filter` | S21/S11 | 用 Q3D C 的模型 vs HFSS：ω_p、κ_p；擬合（ΔL_filter、C_p 比例）→ C_p 擬合值 vs Q3D |
| 3. HFSS `full` | S21/S11 | C_p 固定為第 2 步的值；擬合 ΔL_res、C_j 比例 → C_j 擬合值 vs Q3D，J（Q3D）vs J（擬合） |
| 4.（選做）掃描 | C_j gap 3 個值、C_p plate 3 個值 | J ∝ C_j、κ_p ∝ C_p² 是否成立 |

每次 run 只輸出：電容（Q3D 與擬合值）、HFSS S21/S11、一張圖（HFSS 原始資料 / Q3D C 的模型 / 擬合後的模型）。
存放：`simulation/results/<YYYY-MM-DD>/<HHMM>_<tag>/`，不覆蓋。

## 目前的預測（Design5 整條 channel Q3D：C_p 45.91 fF、C_j 3.43 fF；第 1 步會換成本幾何的值）

| | ω_p (GHz) | ω_r (GHz) | J (MHz) | κ_p (MHz) |
|---|---|---|---|---|
| filter only | 6.7202 | — | — | 37.4 |
| full | 6.7190 | 6.7000 | 19.3 | 37.3 |
