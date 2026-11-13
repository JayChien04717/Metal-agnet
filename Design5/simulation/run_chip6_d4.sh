#!/bin/bash
# Design5 chip (Design4-style channels): one Q3D (8 passes) + one eigenmode per channel (6 passes, 10 um mesh).
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
D4=/c/Users/zx123/Desktop/Metal-agnet-main/Design4/simulation
OUT=/c/Users/zx123/Desktop/Metal-agnet-main/Design5/simulation/results/2026-11-12/chip6_d4
LJ=/c/Users/zx123/Desktop/Metal-agnet-main/Design5/simulation/results/chip6_d4_lengths.json
mkdir -p $OUT
cd $D4
$PY q3d/q3d_hammer_gap.py --spacing 0.4 --gaps 42 --out $OUT/q3d > $OUT/q3d.log 2>&1
for ch in 1 2 3 4 5 6; do
  read LF LR F < <($PY -c "import json;c=[c for c in json.load(open('$LJ'))['channels'] if c['ch']==$ch][0];print(c['L_filter'],c['L_res'],c['f'])")
  MINF=$($PY -c "print(round($F-0.2,3))")
  D4_LENGTHS="$LF,$LR" $PY eigenmode/hfss_hammer_straps.py --cases straps --trims 0 --spacing 0.4 --gap 42 \
     --minf $MINF --modes 2 --out $OUT/eig_ch$ch > $OUT/eig_ch$ch.log 2>&1
  grep "per-point\|modes " $OUT/eig_ch$ch.log | tail -2
done
