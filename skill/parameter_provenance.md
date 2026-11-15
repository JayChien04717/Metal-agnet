# 參數來源與版本規範

適用於新 Design、新 run 及舊結果整理。舊原始資料不覆寫；可用 sidecar 補充來源，不補造缺失證據。

## 角色與狀態

- role：target（目標）、input（輸入）、derived（推導／擬合）、result（計算結果）。下游 input 必須保留上游來源。
- status：draft、recorded（保存不代表驗證）、validated（有驗證證據）、unknown、not_applicable。
- applicability：current、historical、needs_revalidation、unknown。狀態與幾何適用性分開。

## 來源邊界

| 參數 | 來源與要求 |
|---|---|
| EJ、EL | 使用者設定，或量子模型依明確目標／約束求解；記錄固定或可調 |
| lumped Lj／Larray | EM 輸入；保留與 EJ／EL 的換算方向 |
| 電容矩陣 | Q3D 原始結果；記錄 net 順序、參考地、單位與消元 |
| EC、charge coupling | 矩陣推導或 EM 電路擬合；不能混稱直接 EM 輸出 |
| 線性模態、S／Z | HFSS；記錄邊界、ports、網格、頻段及收斂 |
| qubit f01、alpha、ZZ | 量子模型；記錄 EJ／EL／EC、偏壓與截斷 |
| fidelity、leakage | 動力學計算；記錄脈衝、步長、演算法、耗散與定義 |

HFSS 線性頻率不是 qubit 非線性 f01。HFSS 擬合 EC 後，以 scqubits 反求 EJ，記作「scqubits target solve，依賴 HFSS-fitted EC」，不能寫 HFSS 得到 EJ。

使用者固定 EJ／EL 時禁止自行 retarget。歷史授權不明標 unknown。transmon 模型若沒有獨立 EL，記 not_applicable，不以 0 冒充結果。

## 每筆參數

必填 value、unit、role、source_kind、source_ref（確切檔案及欄位／函式）、method、depends_on、constraints、status、applicability、geometry_revision、run_id。未知值可為 null，但需說明原因。

E/h 的 GHz 必須明確標示；f 與 omega、kappa 與 kappa/2pi 不混用。coupling 註明 Hamiltonian 係數定義，不把 charge coefficient 當 exchange coupling。Phi/Phi0 為無因次。

## 每次 run

唯一 run ID 下保存 manifest.json、inputs/、raw/、analysis/。manifest 記錄時間與時區、命令、工作目錄、環境／版本、幾何版本、commit／dirty 狀態、輸入 hash、上游 run、執行狀態及驗證證據。

單有 commit 不足以代表 dirty 工作區，需保存實際輸入快照。重新分析另建 analysis ID 並引用原 run，不覆寫已採用資料。模板未知欄位保持 null；正式 run 執行前填齊，無法補證的舊 run 明確列出缺口。

## 採用結果

design.json 指定 adopted_results，不以最新檔名挑選。README 顯示四捨五入摘要，精確數值回連來源。

幾何、junction 模型或邊界變更時標 needs_revalidation；只有重跑或具體適用性檢查證據才能恢復 current。舊報告的無影響聲明不單獨構成驗證證據。


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
