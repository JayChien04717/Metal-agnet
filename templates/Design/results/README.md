# 模擬結果

有結果時才建立 general／qubit／resonator／flux 分類。每次 run 保存 run.json、code/、raw/，以及需要的 summary.json、圖與 log。使用共用 code/general/run_record.py 建立唯一目錄；入口須明確將輸出寫入該目錄。

重新分析另建 run 並填 upstream；不覆寫原始資料。design.json 明確選定採用結果。
