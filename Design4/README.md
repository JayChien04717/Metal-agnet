# Design4 — Purcell filter + readout resonator

## 可直接匯入的 PurcellResonator

`PurcellResonator` 封裝 10/6 µm CPW 的 Purcell filter、readout resonator 與 feedline。
兩側 hammer 電極和引線已移除，保留開路短支線和原有的 HFSS 跨地線配置。
由 `feed_in`、`feed_out` 接到宿主 layout 的走線或端口。

```python
import sys
from pathlib import Path

from qiskit_metal import designs

# 從 SQL HFSS 資料夾執行，或使用 Design4/pattern design 的絕對路徑。
sys.path.insert(0, str(Path("Design4/pattern design").resolve()))
from purcell_resonator import PurcellResonator

design = designs.DesignPlanar()
design.chips.main.size.size_z = "-500um"
pair = PurcellResonator(design, "readout_channel", options={
    "cpw_spacing": "400um",
    "resonator_trim": "0um",
})
```

- `cpw_spacing` 是兩條平行 CPW 的中心距。預設為 400 µm，小於 250 µm 會顯示警告。無 hammer 版的 HFSS eigenmode 最低實際掃點分裂為 250 µm 對應 2.56 MHz、300 µm 對應 1.34 MHz、350 µm 對應 4.10 MHz。400 µm 參考版為 5.06 MHz。警告門檻採這組離散掃點中符合約 5 MHz 目標的最短間距。掃描資料在 `simulation/results/2026-11-11/no_hammer_spacing/`。
- `filter_length` 和 `resonator_length` 控制兩條走線的名義長度，`resonator_trim` 只移動 resonator 的末端。預設長度來自間距比較用的模型，頻率須依實際宿主晶片和邊界條件驗證。
- `ground_wire_bonds` 預設為 `True`，會標記 HFSS 自動放置跨地線的路徑。晶片尺寸、材料及實體製程跨橋由宿主設計設定。
- HFSS 間距比較結果在 `simulation/results/2026-11-10/purcell_6500_paper_hammer/` 的 `no_hammer_eigen` 和 `no_hammer_400_eigen`。

## Hammer 的地必須接通（2026-11-11）

`stub_f` + hammer pocket + `stub_r` 的蝕刻連成一條，把 filter 與 resonator 之間的地切成上下兩塊。
HFSS 因此多出一條地電流耦合：190 µm、沒有 hammer 電極但 pocket 連續時分裂仍有 122 MHz；
這條路徑不是電容，Q3D 也看不到，擬合會把它全部算進 C_j（10–11 fF，Q3D 為 2.4–2.6 fF）。
Wallraff 在 hammer 兩側 T-junction 放 air bridge（Fig. 2.29 的 C_AB）就是為了這個。

- `purcell_6500.hammer_ground_straps(design)`：每條 stub 中間一條跨地線（跨距 60 µm、線徑 15 µm），落點自動檢查不在蝕刻區。
- `purcell_6500.draw_hammer_straps(renderer, design)`：在 `renderer.render_design(...)` 之後呼叫，畫成 PEC bond wire `g_wb`（與 qiskit-metal auto_wirebonds 相同），Q3D 仍視為地。
  qiskit-metal 的 auto_wirebonds 只在 >400 µm 的直線段放跨地線，所以這兩條 stub 之前一直沒有。
- `build(..., spacing_mm=0.25)` 可直接改間距。示意圖：`fig/purcell_6500_hammer_straps.png`。

**結論：每一段 CPW 都要打跨地線**（hammer stub、T-junction、往 hugger 的 stub、meander 短段都要）。
250 µm、gap 9 µm 的 HFSS eigenmode：只有 qiskit auto bonds 時 J 68.2 MHz；hammer 加 2 條 22.8 MHz；
每段 CPW 都加（62 條，`ground_straps(design, 'cpw')`）20.4 MHz（等效 C_j 3.04 fF，Q3D 2.40 fF）。
腳本 `simulation/eigenmode/hfss_hammer_straps.py --straps cpw`；細節與數據在 `notes/purcell_design4.md`。

## 原 Design6 的驗證記錄


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

**HFSS 流程（Purcell filter 專案的規則）**：先用 eigenmode 切 mesh（first order、線 10 µm / 電容 5 µm 初始 mesh、
Delta F 0.1 %），driven 用 Mesh Link 直接使用那個 mesh，只解 1 個 pass、不再 refine，然後 sweep
（`simulation/driven/hfss6_meshlink.py`）。filter-only 的 HFSS 不需要：沒有 C_in 時 full channel 的 S21 就能同時擬合
ΔL_filter、ΔL_res、C_p、C_j。

每次 run 只輸出：電容（Q3D 與擬合值）、HFSS S21/S11、一張圖（HFSS 原始資料 / Q3D C 的模型 / 擬合後的模型）。
存放：`simulation/results/<YYYY-MM-DD>/<HHMM>_<tag>/`，不覆蓋。

## 目前的預測（Design5 整條 channel Q3D：C_p 45.91 fF、C_j 3.43 fF；第 1 步會換成本幾何的值）

| | ω_p (GHz) | ω_r (GHz) | J (MHz) | κ_p (MHz) |
|---|---|---|---|---|
| filter only | 6.7202 | — | — | 37.4 |
| full | 6.7190 | 6.7000 | 19.3 | 37.3 |


## Design notes

設計紀錄（調參歷史、Q3D / HFSS 數據、程式實作）：[`notes/purcell_design4.md`](notes/purcell_design4.md)。
