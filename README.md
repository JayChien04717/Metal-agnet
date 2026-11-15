# 超導量子晶片 Layout 專案

以 Qiskit Metal 建立版圖，Q3D／HFSS 取得電磁資料，scqubits 建立量子模型，QuTiP 或明確記錄的數值方法分析動力學。實際工具以各 Design 的紀錄為準。

## 開始使用

1. 閱讀 [AI 工作規則](AGENTS.md) 與 [參數來源規範](skill/parameter_provenance.md)。
2. 選擇 Design，確認目前版本、結果適用性與執行環境。
3. 新設計複製 [Design 模板](templates/Design/README.md)，填寫輸入與來源後再執行。

## 內容放置

| 位置 | 責任 |
|---|---|
| `pattern design/` | 共用 QComponent；GDS 工具已移到 code/general/ |
| `code/` | [共用程式與 run 存檔](code/README.md) |
| `skill/` | 必須遵守的流程與驗收要求；[規則入口](skill/skill.md) |
| `note/` | 通用原理、計算方法、技巧與參考資料 |
| `Design*/` | 該設計的輸入、程式、模擬結果、GDS、fig 與歷史 |
| `templates/Design/` | 新 Design 統一模板；未填值不代表可執行 |

skill 寫怎麼做及驗收要求；note 寫原理；設計數值與取捨留在該 Design。README 是摘要，精確結果回連到指定 run。

## 設計索引

[Design1](Design1/README.md) · [Design2](Design2/README.md) · [Design3](Design3/README.md) · [Design3_2](Design3_2/README.md) · [Design3_3](Design3_3/README.md) · [Design4](Design4/README.md) · [Design5](Design5/README.md) · [Design6：新版編排示範](Design6/README.md)

## 工作流程與來源

設計目標／使用者約束 → 幾何與量子輸入 → Q3D／HFSS → 電路推導／擬合 → scqubits → 動力學分析。

EJ、EL 由使用者指定，或由有明確目標與約束的量子模型設計求得；不是 HFSS 的直接輸出。EM lumped L 是輸入。幾何改變後，舊結果仍屬於舊版，不能自動視為新版已驗證。

## 環境與相容

現有 layout／EM 使用 qmetal；Design6 量子分析使用 scqenv。各 run 記錄實際 Python、套件與 solver 版本。本次未建立環境鎖定檔。

所有 Design 已統一為 `layout/`、`simulation/{general,qubit,resonator,flux}/`、`results/{general,qubit,resonator,flux}/`。只有使用中的類別需要建立。混合歷史資料保留在 general/legacy，以各 Design 的 results/index.csv 查詢。

[整理紀錄](note/project_organization.md) · [新 Design 模板](templates/Design/README.md)。既有參數與幾何未調整；design.json 是採用版本索引，不是自動控制所有舊腳本的參數檔。

commit／push 前必須依 AGENTS.md 完成最終版本與參數確認。README 只呈現確認版本的正式資料，舊數據留在歷史紀錄。
