# Design6 — 個別 Purcell 濾波器(Heinsoo 型),依 Cai et al. 2024 做電容結構與 CPW 選擇

參考:Xiao Cai et al., "Optimize Purcell filter design for reducing influence of fabrication variation", Chin. Phys. B 33, 068501 (2024), DOI 10.1088/1674-1056/ad3345。
另:Heinsoo et al., PR Applied 10, 034040 (2018)(濾波器拓樸與公式)。

## 論文重點(摘自全文)
- 拓樸 = Heinsoo:兩個同頻 λ/4 分佈式共振腔(Purcell 濾波器 + readout resonator),電容 C_in(濾波器輸出端,讓反射訊號多數到 P2)、C_p(feedline–濾波器耦合,決定濾波器線寬)、C_j(濾波器–readout resonator,決定 J)。
- 設計值例:C_in = C_p = 40 fF、C_j = 2 fF;f_r = f_p ≈ 6 GHz;κ_p/2π ≈ 45 MHz、J/2π ≈ 8.8 MHz(S21 兩個 dip 分裂)。
- 容差分析:C_p ±4 fF(±11 %)就會破壞 S21 對稱、移動濾波器頻率;C_j ±1 fF(±50 %)只改兩個 dip 間距。
- 電容結構:**T 型包覆電容(開口朝向濾波器)最不敏感**;優於 C 型(開口朝向傳輸線)與 interdigital(IDC 邊緣蝕刻靈敏度約 4.5 fF/µm,特徵尺寸 L=12 µm 時降到 2.1 fF/µm)。
- CPW 選擇:50 Ω 的 10/5.6、13.2/7.5、**16/9 µm**;越大越不敏感;16/9 + T 型(L=12 µm)在 ±0.5 µm 邊緣偏差下,有效耦合 κ_eff 波動只有 ~10 %,r–p 失諧 < 10 MHz。
- 矽基板 ε ≈ 11.5。

## 對本專案的意義
- 目前 Design4/5 的 IDC 對製程敏感;若要實際製作,Design6 用 T 型包覆電容 + 16/9 µm CPW 重做 C_in / C_p / C_j。
- 也是「每個 resonator 一個濾波器」→ 避免共用濾波器的 resonator 之間互相拉動(Design4 看到的通帶畸變)。

## 待辦(尚未開始)
1. 確認 Design5(Jeffrey 型)結果後,決定 Design6 要用哪個 Cj/Cp/Cin(40/40/2 fF 的 Heinsoo 值,頻率改 6.6/6.7/6.8)。
2. 建立 T 型包覆電容元件(qiskit-metal QComponent),16/9 CPW,Q3D 校準電容。
3. 三組(濾波器 + resonator)各自接到同一條 feedline;HFSS driven(15 passes、Interpolating)看 6 個 dip/雙峰。

## 進度(單通道,CPW 10/6,HFSS driven 15 passes、Interpolating 6.3–7.1 GHz)
- 結構已改成論文 Fig. 3(e) 排法:饋線 → C_in(包覆電容,串在線上)→ 節點 → C_p(包覆電容)→ 濾波器(λ/4,往下)與 resonator(λ/4,反向並排),C_j(小電容)跨在兩條平行直線之間,離各自開路端約 1.0 mm。包覆電容形狀依論文 Fig. 3(d) 放大圖量得:厚外框 324×172 µm 挖槽,內有 193×41 µm 的 T 形電極,間隙 ~8 µm,外框離地 ~12 µm。
- Q3D:包覆電容(板 193 µm)= 51.7 fF(目標 40);板縮到 140 µm 後 C_p = 41.5 fF(中心板側對地 14 fF,外框側對地 113 fF)。C_j(2 指、錐 40 µm)的 Q3D 連續三次 HFSS COM 錯誤(舊版 4 指/2 指曾量到 4.39/1.52 fF)。
- HFSS:舊的串聯結構完全平(S21 −15 dB,S11 0 dB)。新結構 c6_v3(板 193、濾波器 4.227 mm、resonator 4.484 mm):6.4005 GHz 有 −42 dB 的深 dip 加 Fano 峰、6.544 GHz 有 0.9 dB 小 dip;c6_v4(板 140、濾波器 −0.125 mm、resonator −0.105 mm):只剩 6.3455 GHz 的 −62 dB dip。
- **矛盾**:ABCD(點式模型)預測 v3 的 dip 在 6.60/6.69 GHz、v4 在 6.84/6.89 GHz(縮短應往高頻移 +240 MHz),HFSS 的 dip 反而從 6.40 往下移到 6.35 GHz。目前不知道這個 dip 是什麼(可能是 C_in/C_p 與饋線的寄生模態,不是濾波器或 resonator)。
- 待辦:分段 HFSS 驗證——先只放 C_in + C_p + 濾波器(無 resonator、無 C_j),量濾波器 dip 的實際頻率與線寬,再加 resonator;修好 C_j 的 Q3D(或改用論文的橫向小電容形狀)。

## 濾波器單獨驗證(f6_v2,C_in + C_p + 濾波器,無 resonator / C_j,6.3–7.1 GHz Interpolating 15 passes,20 min)
- 單一 dip:**6.653 GHz,−58.5 dB**(濾波器共振,目標 6.7 GHz,低 0.7 %);基線最高 −14.8 dB(6.70 GHz 附近,C_in 串聯造成)。
- 因此濾波器本身在預期頻率附近;之前完整通道的 6.35 GHz dip 不是濾波器單獨造成的 → 下一步加 resonator(先不放 C_j 看其單獨 dip),再加 C_j。
- DRC(skill/chip_drc.md):resonator 蛇形改為 spacing 350 µm / YM 1.2 mm(最短段 350 µm,原 67 µm 違反 B2);EnclosingCap 板與外框尖角倒圓 8 µm、pocket 倒圓(C1);plot_design 畫出 99 µm 圓角。T 形接點仍為直角。電容改了,Q3D 需重量。

## 完整通道(c6_v5,DRC 修正版,32 min)
- dip 6.544 GHz(−61.0 dB,深)與 6.681 GHz(−16.0 dB,prominence 0.7 dB,淺);基線最高 −14.0 dB。
- 對照:濾波器單獨 6.653 GHz;舊版完整通道 6.3455 GHz。加入 C_j + jf 分支 + resonator 後濾波器 dip 下移約 110 MHz;6.681 GHz 的淺 dip 很可能是 resonator(弱耦合,未證實)。

## 紀錄與比較(每次模擬都存 S11/S21 與對應電容)
- `simulation/run_record.py`:`hfss_chan6.py` 每次跑完自動寫 `results/<tag>_dB.csv`(S11/S21 dB)與 `<tag>_meta.json`(Q3D 最新的 C_in/C_p/C_j、板尺寸、C_j、蛇形參數、長度、備註;備註用環境變數 `RUN_NOTE`)。`q3d_caps6.py` 每次也存一份帶時間的 `q3d6_<cap>_cmat_fF_<MMDD_HHMM>.csv`。
- `simulation/compare_runs.py [tag...]`:疊圖 `results/compare_S.png` + 表格 `results/runs_summary.csv`(dip 位置、各電容)。
- 舊 run(c6_v1/v3/v4/v5、f6_v1/v2)已用 `backfill_records.py` 補記;電容為當時量到的值,沒量到的留空(C_in 圓角版尚未重量)。f6_v1 不可靠(誤啟動的舊程式完整通道)。
- C_j Q3D:taper 40 µm 的版面一律 COM 錯誤;taper 100 µm 可算(overlap 12 µm → 1.52 fF),overlap 40 µm + taper 100 也失敗。

## 清理與改版(依「只保留論文設計、C/T 型、hammer」)
- 已刪除:Design5(Jeffrey 型,整個資料夾)、Design6 舊版 run(c6_v1/v3/v4/v5、f6_v1/v2 的 S 參數、meta、圖、log)、`design6_layout_v1_serial.py`、舊對照圖與舊 GDS、補記腳本。保留 Design1–4(Design6 要 import 其中 ftf_layout、purcell_layout、interdigital_cap、tap_junction、res_end)。
- 論文 Fig. 3(Cai 2024,存為 `fig/paper_cai2024_page_fig3.png`):T 型 = 包覆電容開口朝向濾波器;C 型 = 開口朝向傳輸線(圖上 (c)/(d) 標籤與圖說相反,以內文為準);C_j = 「兩個矩形」互相嵌合的小電容(我們叫 hammer 型,實作為 2 指、無 taper 的 InterdigitalCap),兩條垂直線距離約 130–150 µm 量級。
- 版面開關:`CIN_TYPE`、`CP_TYPE` = 'T' | 'C'(C 型 = 同一元件反向接,框接濾波器側);兩種都能建出、轉角與線段檢查通過(`fig/design6_layout_T.png`、`_C.png`)。
- C_j:無 taper,overlap 25 µm → 2.00 fF(24 µm 1.95、40 µm 2.71),對地各 15 fF(taper 100 版為 33 fF);兩線距 0.30 mm。C_in 40.57、C_p 40.61 fF。
- 完整通道 HFSS(T 型)c6_T1 跑中。

## hammer 電容(`pattern design/hammer_cap.py`,類別 `HammerCap`)
- 依論文 Fig. 3(e) 的 C_j 放大圖重寫成獨立元件:CPW 把手 → 矩形「錘頭」(14 × 30 µm)→ 各一根手指(12 µm 寬)交錯越過對方,兩頭相距 17 µm,蝕刻口袋 6 µm,頭部圓角 2 µm。取代之前拿 InterdigitalCap 湊的細梳(那個頭只有 18 µm 高,和論文不像)。
- Q3D:C_j = 1.97 fF(目標 2.0),對地每側 13.3 fF;overlap 5 µm。兩線距 0.30 mm。
- 圖:`fig/zoom_cj.png`(C_j 放大),論文放大圖 `fig/_cj_inset-4.png`。

## hammer 修正版(讀論文 Fig. 3e 的 C_j 原生解析度圖後重做)
- 論文 C_j 放大圖(CPW 16/9,50 µm 比例尺):橋(CPW 把手)兩端各接一個**矩形板(錘頭)**,兩板**互相平行**、垂直於橋,中間一條窄縫;板 ≈ 25 × 65 µm、縫 ≈ 15 µm、蝕刻口袋比板外擴 ≈ 16 µm(≈ 96 × 97 µm);兩條長線的間距約 100 µm(邊到邊)。之前的版本(交錯手指、細梳)都不是這個形狀。
- `hammer_cap.py` 重寫成這個形狀:`head_width 25 / head_height 65 / head_gap / lead_length 10 / ground_gap 16`,頭部圓角 2 µm。
- Q3D(10/6 CPW):縫 15 µm → 3.06 fF;25 µm → 2.07;35 µm → 1.48;**縫 26 µm → 2.00 fF**(對地每側 13.3 fF)。兩條垂直線距改為 0.20 mm。
- Q3D COM 錯誤的原因:導線路徑起點伸進板內(重疊),改成貼著板邊(1 µm)就能算;與 HFSS 關掉重開無關(重開後 C_p 仍是 40.61 fF)。
- 完整通道 HFSS:c6_T2(T 型 C_in/C_p + hammer C_j)。

## c6_T2 結果(T 型 C_in/C_p + 平行板 hammer C_j,10/6 CPW)
- 兩個 dip:**6.585 GHz(−49 dB,濾波器)**、**6.731 GHz(−16 dB,淺,應為 resonator)**;基線最高 −14.7 dB。濾波器單獨 6.653 GHz → 加入 C_j/分支/resonator 後下移 68 MHz。
- 兩者差 146 MHz,論文要求 ω_r = ω_p。校正:濾波器縮短 0.077 mm(−0.125 → −0.202),resonator 加長 0.02 mm(−0.105 → −0.085)→ c6_T3。

## 理論曲線 vs HFSS(`simulation/theory_abcd.py`)與結果資料夾
- 每次模擬結束後自動把該次檔案放進 `results/<YYYYmmdd_HHMM>_<tag>/`:`<tag>_S.csv`(複數)、`_dB.csv`(HFSS S11/S21 dB)、`_meta.json`(電容、板尺寸、長度、備註)、`_theory_dB.csv`、`_theory_params.json`、`_theory_vs_hfss.png`(HFSS 實線、理論虛線)、`_q3d_<cin|cp|cj>_cmat_fF.csv`(當時的 Q3D 矩陣)。不同次模擬用時間戳分開。
- 理論模型(ABCD/導納):P1 – 線 – C_in(π:Q3D 串聯 C + 兩電極對地 C)– 線 – 節點 n0(並聯分支)– 線 – P2;分支 = 線 – C_p(π)– 濾波器線(C_p 端開路、遠端短路)– 離 C_p 1.02 mm 處 C_j(π)– resonator 線(開路端 d_r、短路端 L_res)。Z0 = 50.9 Ω,v = 1.2016e8 m/s,無損。
- `compare_runs.py`:所有 run 疊圖(HFSS 實線 + 理論虛線)+ `results/runs_summary.csv`。
- c6_T2:理論 S21 基線 −18 dB(HFSS −17 dB,接近);理論 dip 6.654 GHz(濾波器)與 6.740 GHz(resonator,窄);HFSS 6.585 與 6.731 GHz。理論濾波器 dip 與 HFSS 濾波器單獨(6.653)幾乎相同,HFSS 加上 C_j/分支後多下移 ~70 MHz,理想模型沒有這個效應(橋/接點的分佈效應,待查)。

## 收斂資料集 + 理論擬合(2026-11-08)
- C 槽滿(0.1 GB)造成 HFSS COM 錯誤;刪除 Temp 內三個舊 hfss_*.pjt(86 GB)。建議 HFSS Temp 目錄改到 D:。`hfss_chan6.py` 修正語法錯誤,新增 `MAXP` 環境變數。
- c6_T3h(T3 幾何,solve 6.62 GHz,6.45–6.80 Interpolating 0.25 MHz,18 passes,137 min):單一 dip **6.6182 GHz**(−64 dB)。ΔS 未收斂(0.079;在 dip 上求解時 ΔS 在 0.08–0.2 間跳動),但 dip 位置 T3/T3g/T3h = 6.623/6.619/6.618,mesh 相依約 5 MHz。
- c6_T3i(加 CAP_MESH 12 µm):2 h 40 min 後 COM 錯誤(非磁碟,可能記憶體)。
- e6_T3b eigen(15 passes,ΔF 0.04 %):**6.3419、6.6670 GHz**。
- 擬合(`fit_theory.py c6_T3h --eig 6.3419,6.6670`,方法見 note/purcell_filter_method.md):dLf +0.070、dLr +0.154 mm、s_in 0.98、s_p 0.82、s_j 0.80、kv 0.978;rms |ΔS21| 0.081 → 0.017;兩個 eigen 頻率與 driven dip 都吻合(6.6189 vs 6.6182)。
- 結論:**6.619 GHz = 濾波器,resonator 在 6.342 GHz**(dLr +20 µm 時移動的是 6.342),失諧 277 MHz,所以 driven 只看到一個 dip(resonator 經濾波器耦合到線上的線寬只有 kHz 級)。
- 模型給出:κ_p/2π ≈ 32 MHz(目標 45)、J/2π ≈ 11.9 MHz(目標 8.8;J 在大失諧下無法由資料決定,主要來自 Q3D C_j)。
- 下一步:resonator 縮短約 0.19 mm(L_RES_TOTAL 4.3986 → ~4.21)讓 ω_r ≈ ω_p,再跑 driven 6.50–6.75 GHz,屆時兩個 dip 分開 ~2J,κ_p 與 J 可直接擬合。
- 最小 step / pass 整理(詳見 skill/skill.md §10):最小 step = c6_T3h 的 0.25 MHz(Interpolating);driven 最少 15 passes 即可(15→18 只差 0.8 MHz、多 ~100 min),eigen 最少 14 passes(與 15 差 < 1 MHz)。

## 做成 paper 的 S21 形狀(c6_T5 版面,2026-11-08,尚未模擬)
- 方法(`simulation/design_target.py`,用 c6_T3h + e6_T3b 擬合後的模型):
  κ_p = 濾波器負載極點的線寬(電路複數零點,C_j → 0);目前幾何只有 22.5 MHz(之前 notch 寬 32 MHz 被傾斜背景誇大)。
  paper 的兩個 dip = S21 傳輸零點 = n0 接地時濾波器 + resonator 的無損模態;改 resonator 長度時最小間距 = 2J(與 κ_p 無關)。因 κ_p/4 > J,帶損耗的混成模態不分裂,不能用來定 J。
  依序:C_p → κ_p 45 MHz;濾波器長度 → 極點頻率;resonator 長度 → 兩零點反交叉;C_j → 最小間距 17.6 MHz(J 8.8);外層迴圈讓兩零點中心 = 6.70 GHz。
- 不對稱來源:C_in 到 n0 的線(C_in 當反射面形成小腔,Fano)。l_feed_c 1.25 → 0.25 mm 時 peak 0.185 → 0.161,接近右基線 0.155(`results/asym_check.png`)。
- 結果(LFC 0.25,`results/design_target_lfc0.25.png/json`):零點 6.691 / 6.709 GHz、κ_p 45.0、J 8.8 MHz、peak ≈ 基線。原位 C_in 需要 C_p 60 fF 且 peak 0.19 > 基線,不採用。
- 版面變更:C_in 移到 CIN_X −0.3025(C_in→n0 線 20 µm);C_p 改用 ENC_P 板寬 140 → 175 µm(目標 Q3D 47.3 fF,用斜率 0.19 fF/µm 估,**需 Q3D 重量**);C_j 間隙 26 → 35 µm(Q3D 1.48 fF,目標 1.49);濾波器 offset −0.202 → −0.2975 mm;resonator −0.276 → −0.329 mm;XM 0.9 → 0.87 讓蛇形腳 355 µm(≥ 350 規則)。
- 圖:`fig/design6_layout.png`、`fig/design6_T5_zoom_cin_cp.png`、`fig/design6_T5_zoom_cj.png`。
- 下一步:Q3D 重量 C_p(175 µm)與 C_in(新位置旁有 C_p),再跑 HFSS c6_T5(solve 6.70,6.60–6.80 GHz,0.25 MHz,MAXP 15)。

## c6_T5 HFSS 結果(2026-11-08)— 沒有做出 paper 形狀
- Q3D:C_p 47.28 fF(175 µm,命中)、C_j 1.48 fF(35 µm)、C_in 40.57 fF;C_p 對地 121.6/15.8 fF(模型假設 ~131/17)。
- c6_T5(6.60–6.80 GHz Interpolating 0.25 MHz,15 passes,32 min):整段平坦 −17 dB → **漏掉了 dip**(Interpolating 收斂但沒取樣到共振)。
- c6_T5w(5.8–7.6 GHz,1 MHz,12 passes,20 min):兩個寬 dip **6.276 GHz(FWHM 64 MHz)、6.754 GHz(32 MHz)**。
- f6_T5w(只有 C_in + C_p + 濾波器,寬掃):單一 dip **6.836 GHz** → 第二個寬 dip 需要 C_j / resonator 分支。
- 幾何比對(HFSS 物件 T3h vs T5w):只有預期的改變,C_j 兩板相距 35 µm,沒有短路。
- c6_T5r(resonator +0.1 mm,DLR_EXTRA):兩個 dip **都往上** 移約 35 MHz(6.311 / 6.791)。加長 resonator 應往下 → 兩個 dip 都不是 resonator;12 passes 寬掃的 mesh 誤差可能有 ~35 MHz。
- 結論:電路模型(只用 C_j 點耦合 + 擬合長度修正)無法外推到 c6_T5;分支(jf/jr 接點、C_j 對地、兩條平行線)在 HFSS 裡的效應比模型大很多,且濾波器單獨就比預期高 ~120 MHz。

## 混合模型(EM block + 電路)與照 paper 的 c6_P1(2026-11-08)
- 直接擬合 c6_T5w/T5r 失敗:模型只能產生一個寬 dip(結構上缺分支的物理),C_j 放寬到 100× 也一樣。
- filter-only HFSS(f6_T5w 6.836)與**原始 Q3D 模型**(6.851)只差 15 MHz → 濾波器部分模型正確;c6_T5 的偏差來自用了 c6_T3h 擬合的錯誤修正係數(s_p 0.82、s_j 0.80、kv 0.978)。
- 混合模型:`simulation/hfss_block6.py`(C_j 分支 3-port:濾波器 C_p pin b → jf、resonator 開路端 → jr、C_j)+ `theory_abcd.load_block/block_Yin` + `hybrid_check.py`。
  - blk_T5:resonator 開路端被接地(ResEnd 只有 pin,沒列在 open_terminations 就會短路)→ 修正後 blk_T5b(20 passes,ΔS 0.0028)。
  - blk_T5b:S33 與理想差 ~4°(OK);**|S31| ≈ 0.017 = 理想 C_j 的 6 倍**,相位差 ~70°,原因未確認(地平面太小 / lumped port / 真的耦合)。
  - 混合模型:濾波器 6.648、resonator 6.924 GHz;resonator +0.1 mm → 6.780(−144 MHz,合理)。HFSS 12-pass 寬掃(T5w/T5r)兩個 dip 一起往上 35 MHz,不可靠。
- **c6_P1(照 paper 的電容)**:C_in = C_p = 40 fF(板 140 µm,Q3D 40.6)、C_j = 2 fF(hammer 縫 26 µm,Q3D 2.00)、C_in 緊鄰 n0(與 paper Fig. 3e 相同)、兩線中心距 0.2 mm(paper 約 0.2 mm)、C_j 離開路端 1.02 mm。
  長度用原始 Q3D 模型:L_FILTER offset −0.1533(4.0734 mm)、L_RES_TOTAL offset −0.0807(4.4029 mm)→ 零點 6.685 / 6.715 GHz,κ_p 48.5 MHz(paper 45.5)、J 14.9 MHz(paper 8.8,頻率較高 + 幾何不同)。圖 `simulation/results/paper_like_model.png`。
- Q3D 檔 q3d6_cp/cj_cmat_fF.csv 已還原成 140 µm / 26 µm 的值(從 c6_T3 資料夾)。
- 未確定:混合模型 block 的 6× 耦合;尚未跑 c6_P1 的 HFSS。
- **c6_P1 HFSS**(Fast 5.8–7.6 GHz,16 passes 收斂 ΔS 0.002/0.0029,40 min):兩個寬 dip **6.172 GHz(FWHM 64 MHz)、6.582 GHz(25 MHz)**,相距 410 MHz;原始模型預測 6.685/6.715。與 c6_T5(6.276/6.754)同樣的「兩個相距很遠的寬 dip」→ 這是真的(已收斂),不是 mesh 雜訊。
  混合模型(用 blk_T5b,縫 35 的 block 近似)給 6.482/6.577,第二個 dip 對上 6.582,第一個差 310 MHz。圖 `simulation/results/c6_P1_vs_models.png`。
  結論:分支區有遠大於 C_j(2 fF)的耦合,Q3D 單獨量 hammer 抓不到。下一步:拿掉 C_j(hammer + stub)跑全通道,看兩個寬 dip 是否還在。

## Wallraff parameters (2026-11-09)

* Reference screenshots and the transcribed Table 2.4: `SQL HFSS/note/wallraff reference/` (README.md).
* `simulation/wallraff_model.py`: their circuit (eqs. 2.122–2.136) with direct pole extraction of ω_p, ω_r, J, κ_p. It reproduces their J^Sim (14–15 vs 15–16 MHz) with one v = 1.147e8 m/s.
* **Main point:** C_pr ≈ 2.5 fF gives J ≈ 15 MHz in their design and in our model. The 410 MHz splitting of c6_P1 in HFSS therefore cannot come from the 2 fF C_j.
* **W1 design (circuit only):**
  * Q35 capacitances with Design6 CPW;
  * filter 3.629 mm (tap 1.2 mm), resonator 4.458 mm (tap 1.3 mm from the open end);
  * J 16.9 MHz, κ_p 37.6 MHz, dips 6.690 / 6.728 GHz.
  * Plot: `simulation/results/wallraff_model.png`.
* **To build W1 in the layout:**
  * flip the hugger (frame on the filter side) and set C_pf to about 45 fF;
  * hammer about 2.8 fF;
  * Q3D both;
  * then HFSS (show the layout first).

### W1 layout (2026-11-09)

* Built with `D6_VARIANT=W1` in `pattern design/design6_layout.py`. The default is still c6_P1.
* Output: `gds/design6_W1.gds`, `fig/design6_W1_layout.png`, and zooms in `simulation/results/W1_zoom_cp.png` and `W1_zoom_cp_cj.png`.
* **Changes from P1:**
  * CP_TYPE 'C' (hugger frame toward the filter, plate toward n0);
  * C_p plate 165 µm (estimated 45.5 fF);
  * C_j gap 19 µm (estimated 2.77 fF);
  * D_CJ 1.20 mm;
  * OPEN_DY 0.87 mm (resonator tap 1.30 mm from the open end);
  * XMF / YMF 0.5 / 0.5;
  * L_FILTER 3.6289 mm, L_RES_TOTAL 4.4578 mm.
* **Checked from the pins:** filter 3.6289 mm, tap 1.200 mm; resonator 4.4578 mm, tap 1.300 mm.
* **Next steps:**
  1. Q3D of cp and cj (W1 geometry).
  2. Re-solve the lengths with the Q3D π networks in theory_abcd. Our C_j ground C is about 13 fF (Wallraff 2.4), and the frame-to-ground C is about 110 fF (Wallraff 93).
  3. HFSS.

### W1 Q3D and re-solved lengths (2026-11-09)

* **Q3D (W1 geometry, `D6_VARIANT=W1 python q3d_caps6.py cp cj`):**
  * C_p = 45.41 fF; frame 119.3 fF to ground on the filter side, plate 15.4 fF on the n0 side;
  * C_j = 2.59 fF (19 µm gap), 13.0 fF to ground on each side.
  * `results/q3d6_{cp,cj}_cmat_fF.csv` are now the W1 values. The P1 values are in `q3d6_{cp,cj}_cmat_fF_P1.csv`.
* **`simulation/solve_w1.py`:** raw Q3D π networks in the theory_abcd model; ω_r / ω_p from the branch zero and the poles.
  * As drawn (Wallraff-cap lengths): ω_p 6.469, ω_r 6.614 GHz. The filter is pulled down by the 119 fF frame.
  * Solved: filter 3.4542 mm (−0.1747), resonator 4.4070 mm (−0.0508).
  * Result: ω_r 6.700, ω_p 6.719 GHz, **J 14.5 MHz, κ_p 36.6 MHz** (Wallraff Q35 sim: 15.0 / 41.0).
  * S21 zeros 6.6917 / 6.7256 GHz.
* The full-channel baseline is about −16 dB. This is the insertion loss of C_in = 40 fF in series (Wallraff eq. 2.121: 44 fF → 14 dB).
* The lengths are now set in the W1 block of design6_layout.py. GDS and figure have been regenerated.
* **Next:** HFSS full channel (Fast, 6.4–7.0 GHz or wider, about 16 passes).
