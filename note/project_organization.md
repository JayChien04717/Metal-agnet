# 專案內容編排

程式按 general、qubit、resonator、flux 分類；layout 與 QComponent 獨立。原始 Design1、Design2、Design3、Design3_2、Design3_3、Design4、Design5、Design6 全部已套用。

共用工具移到根目錄 code/；仍含設計常數與跨設計依賴的程式留在所屬 Design，不以複製或合併掩蓋差異。Design3_2 的製程元件與其他設計共用來源有差異，保留其本地版本。

舊散列檔案按主題放入各類別的 legacy；真正混合批次與彙整保留 general/legacy。已有完整 manifest 的獨立 run 保持內部結構不變。results/index.csv 記錄每個原路徑、新位置與內容 hash，原 run 的 manifest 與 source snapshots 不改寫。

刪除無現役 Python 引用的舊 reference 腳本與已完成的固定日期 batch queue；保留 reusable solver 與原始結果。具體清單見 [搬移與清理紀錄](../code/general/migration_record.json)。20261115_latest 是原有交付快照，本次未更動。

所有正式參數仍需在 commit／push 前依使用者確認；B2 的候選確認保留。目錄搬移不等於重新驗證模擬。

驗證結果：[migration_validation.json](../code/general/migration_validation.json)。8 個設計的搬移前後 QGeometry 一致；746 份搬移資料中 745 份 hash 相同，另 1 份 Design2 manifest 在整理期間由執行中的工作更新，已保留恢復分析前後紀錄。
