# 專案 AI 工作規則

## 任務入口

- 先讀根目錄 README、目標 Design 的 README 與 design.json（若有）。
- 參數／結果工作：讀 skill/parameter_provenance.md。
- QComponent／layout：讀 skill/skill.md、skill/qcomponents.md、skill/chip_drc.md；launch pad 加讀 skill/launch_pad.md。
- GDS／fig：加讀 skill/gds_export.md、skill/layout_figures.md。
- 電容／互感：分別加讀 skill/cap_matrix_analysis.md、skill/mutual_inductance_analysis.md。
- 新 Design 複製 templates/Design/ 並填寫，不把模板當成完成設計。

## 編排

- 共用 QComponent 放 pattern design/；新 Design 組裝放 layout/，專用元件放 layout/components/。
- skill 放操作規則；note 放通用原理；設計數值、例外與歷史留在該 Design。
- 所有 Design 已按用途整理；跨設計匯入與歷史結果路徑由 code/general/project_paths.py 和搬移索引銜接。新增程式直接使用新位置。
- 保留其他未提交變更，不覆寫、不回復。

## 參數與證據

- 使用者限定的 EJ／EL 不得自動 retarget；無法滿足目標時記錄差距。
- EJ／EL 是使用者或量子模型設計的輸入；HFSS lumped L 也是輸入。不得把它們宣稱為 HFSS 獨立求得的量。
- target、input、derived、result 分開；按參數記錄來源、單位、假設與驗證狀態。
- 用 design.json 指定採用的 run／檔案，不靠最新時間選正式結果。
- 缺失的來源、授權、commit、hash 或收斂證據標 unknown，不補造。
- 幾何變更後標記舊結果 historical 或 needs_revalidation。

## 執行與交付

- layout／EM 使用該設計指定的 qmetal；量子分析使用指定的 scqenv 或其他已記錄環境。
- 保留既有 solver 預設。pass 數是執行設定，不等於收斂證據；記錄實際檢查或缺口。
- 新 run 保存不可覆寫的輸入快照、原始資料、分析、環境與狀態；重新分析另存並引用父 run。
- simulation／fab 幾何分別標記；GDS 與 fig 必須對應同一版本。
- README 放目前摘要、來源、重現入口與限制；完整歷史放 notes。
- 文件變更檢查連結與 JSON；未重跑模擬必須明確說明。

本文件與參數來源規範更新舊規則中「所有腳本使用同一環境」及「讀最新參數檔」的約定；其他既有幾何／製程規則繼續適用。


## Commit／push 前的最終版本確認

- commit 或 push 前，先整理可供審閱的版本摘要：最終候選、參數設定、採用的 run、GDS／fig 及未驗證項目，再請使用者確認最終使用版本。
- 使用者已明確確認且版本／參數未改變時沿用確認，不重複詢問；任何後續參數或版本微調都使該次確認失效，提交前須重新確認。
- 版本確認不等同 commit／push 授權；仍按使用者要求執行 Git 操作。
- README 正式參數與結果只對應最後確認的版本。舊版數值放 notes／歷史來源索引；新版缺少結果時寫待驗證，不以舊資料填補。
- design.json 記錄 final_selection（版本、確認來源、確認範圍）及 release_confirmation（完整參數快照是否已確認）。選定候選不代表所有參數、結果或製程版本均已確認。


## 目前目錄與存檔約定（2026-11-15）

以使用者最後選擇的簡化分類為準：共用工具放 code/general 或 code/qubit（有需要再加 resonator／flux）；各 Design 用 layout/、simulation/{general,qubit,resonator,flux}/，結果用同樣分類。不要另加 jobs、support、release 或 sweep 層級。

新 run 用 run.json、code/、raw/ 加摘要／圖／log；設定值與來源寫 run.json，需封存輸入時加 inputs/。重新分析另開 run 並填 upstream。舊混合結果保留 general/legacy，原始記錄不改寫；搬移索引不是歷史來源證據。

部分舊入口仍寫固定檔名，執行前必須接入新 run 的輸出目錄，不能覆寫歷史結果。路徑相容工具不等於所有舊入口已自動接入 run recorder。
