# DesignNN：填寫設計名稱

> 複製整個 templates/Design 到專案根目錄，命名 DesignNN，再替換此文件與 JSON 的佔位內容。此模板沒有可執行版圖或模擬程式；JSON 尚需由設計程式明確載入。空值表示未設定，不能當成 0。

## 1. 設計目的與目前狀態

填寫拓樸、qubit 種類、晶片大小、設計目標。狀態：draft。指定 [design.json](design.json) 的 geometry_revision 與採用結果；未模擬就明確寫未模擬。

## 2. 版圖與產物

完成後加入 fig 全晶片／局部圖及 gds 連結，註明 simulation 或 fab、幾何版本、來源 run。未產生前不建立虛假產物連結。

## 3. 參數摘要與來源

| 參數 | 目標 | 採用值 | 單位 | 來源／狀態 |
|---|---|---|---|---|
| 待填 | 待填 | 未設定 | 待填 | draft |

[幾何輸入](parameters/geometry.json) · [量子輸入](parameters/quantum.json) · [solver 設定](parameters/solver.json) · [來源表](parameters/provenance.json)

EJ／EL 來自使用者設定或明確的量子目標求解；若固定，不可自動 retarget。EC 與 coupling 記錄推導／擬合方法。精確值只從指定的 run 引用，不手動維護多份互相競爭的正式數值。

## 4. 模擬結果與適用範圍

尚無採用結果。完成後每項列出 run ID、原始資料／分析連結、幾何版本、驗證證據與限制。更新 design.json 的 adopted_results；不按檔案時間挑選。

## 5. 重現入口與環境

尚無可執行入口。完成 layout/build.py 與所需 simulation 腳本後，記錄從哪個工作目錄執行、完整命令、Python 環境及工具版本。執行前填齊輸入；每次 run 使用唯一 ID。

## 6. 待完成與目錄

- 填寫設計輸入、來源、固定／可調約束。
- 實作版圖與需要的求解流程；記錄 simulation／fab 的差異。
- 完成 DRC、輸出一致性、EM 收斂與量子截斷檢查。

```text
DesignNN/
  README.md
  design.json
  parameters/       geometry.json, quantum.json, solver.json, provenance.json
  layout/           晶片組裝；components/ 放設計專用 QComponent
  simulation/       general/、qubit/、resonator/、flux/
  results/          按用途分類；每次 run 的 run.json、code、raw 與摘要
  gds/              匯出產物
  fig/              結構圖
  notes/            設計取捨、調參歷史與例外
```

共用規則在專案根目錄 AGENTS.md 與 skill/parameter_provenance.md。歷史設計可以保留路徑，但須在 design.json 明確記錄對應。

## Commit／push 前確認

先整理最終候選、參數快照、run、GDS／fig 與未驗證項目，交使用者確認並記入 design.json。候選選定與完整參數確認分開；後續微調須重新確認。README 只呈現最後確認版本的資料；舊數據移入 notes，新版未驗證就留待驗證。
