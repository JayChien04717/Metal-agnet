# 共用程式

| 類別 | 內容 |
|---|---|
| general | GDS、繪圖、AEDT 連線、路徑相容、run 存檔 |
| qubit | 電容矩陣與量子分析的共用算法 |

resonator／flux 在出現真正共用程式時再建立；不把仍綁定 Design 常數的腳本冒充共用庫。共用 QComponent 仍在根目錄 pattern design/。

`project_paths.py` 與 `path_migration.json` 維持舊資料路徑相容，沒有複製舊目錄或建立捷徑。歷史 manifest 不改寫，可由索引查新位置。

## 新 run

從專案根目錄執行（只準備紀錄，不啟動 solver）：

```powershell
conda run -n qmetal python code/general/run_record.py Design3_3 flux B2_check --entrypoint Design3_3/simulation/flux/recheck_readout_flux.py
```

指令回傳唯一 run 資料夾。將已確認輸入用 --input 傳入，父結果用 --upstream 記錄。solver 入口需明確接收新目錄，不能只建立紀錄卻仍把資料寫回舊結果。

程式內可用 `create_run(..., parameters=實際設定)`，將輸出寫進回傳資料夾，再以 `finish_run(..., status='completed', summary=結果)` 結束；失敗記 failed。收斂未檢查時 validation 保持 unknown。geometry_revision 與 solver_version 必須補實際值，不猜測。

為涵蓋跨 Design 依賴，code 快照保存所有現役 Python 與路徑索引；不包含歷史 results 中的副本。第三方套件記版本，沒有宣稱已封存整個執行環境。
