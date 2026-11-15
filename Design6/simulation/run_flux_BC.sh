#!/bin/bash
# Flux-line variants B (L short: sep 9 um, bar 60 um) and C (T type: sep 9 um, arms 2 x 50 um), local HFSS model,
# then dense sweep and the acceptance metrics (M, 30 ns ring-up, transmon-model crosstalk matrix, square pulse).
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
export AEDT_PROJECT="$(cd "$(dirname "$0")" && pwd -W)/aedt/Design6_flux.aedt"   # Design6.aedt would not open after the 2026-11-14 crash
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
SQ=/c/Users/zx123/miniforge3/envs/scqenv/python.exe
R=results/2026-11-14
declare -A OV
OV[Lshort]='{"coupler": {"flux_bias_line_options": {"make_fbl": true, "fbl_sep": "9um", "fbl_length": "60um", "fbl_x": "-25um", "fbl_type": "L"}}, "qflux": {"fbl_sep": "9um", "fbl_length": "60um", "fbl_type": "L"}}'
OV[T]='{"coupler": {"flux_bias_line_options": {"make_fbl": true, "fbl_sep": "9um", "fbl_length": "100um", "fbl_x": "-25um", "fbl_type": "T"}}, "qflux": {"fbl_sep": "9um", "fbl_length": "100um", "fbl_type": "T"}}'
for V in Lshort T; do
  $PY -I driven/hfss_flux.py --tag loc_$V --overrides "${OV[$V]}" > results/hfss_flux_loc_$V.log 2>&1; echo "$V hfss exit $? $(date +%H:%M)"
  $PY -I driven/flux_resweep.py --tag loc_$V --dense > results/flux_dense_$V.log 2>&1; echo "$V dense exit $? $(date +%H:%M)"
  $SQ -I driven/flux_metrics.py --run $R/flux_loc_$V --tag $V > results/metrics_$V.log 2>&1; echo "$V metrics exit $?"
  $SQ -I driven/flux_xtalk_matrix.py --run $R/flux_loc_$V --tag $V > results/xtalk_$V.log 2>&1; echo "$V xtalk exit $?"
  $SQ -I driven/flux_square_spectro.py --run $R/flux_loc_$V --line Q1 --amp 0.24 --tag $V > results/sq_Q1_$V.log 2>&1; echo "$V square exit $? $(date +%H:%M)"
done
echo ALL DONE
