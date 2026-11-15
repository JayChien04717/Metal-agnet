#!/bin/bash
# Flux-line variant B2 (2026-11-15 layout: L, sep 9 um, bar 45 um, fbl_x 12 um for Q1 / Q2 / TC), local HFSS model,
# dense sweep, acceptance metrics, then the A-D + B2 comparison. Results go to results/2026-11-14/ next to A-D.
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
export AEDT_PROJECT="$(pwd -W)/aedt/Design6_flux.aedt"
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
SQ=/c/Users/zx123/miniforge3/envs/scqenv/python.exe
D=2026-11-14; R=results/$D; V=B2
OV='{"coupler": {"flux_bias_line_options": {"make_fbl": true, "fbl_sep": "9um", "fbl_length": "45um", "fbl_x": "12um", "fbl_type": "L"}}, "qflux": {"fbl_sep": "9um", "fbl_length": "45um", "fbl_x": "12um", "fbl_type": "L"}}'
$PY -I driven/hfss_flux.py --tag loc_$V --date $D --overrides "$OV" > results/hfss_flux_loc_$V.log 2>&1; echo "$V hfss exit $? $(date +%H:%M)"
$PY -I driven/flux_resweep.py --tag loc_$V --date $D --dense > results/flux_dense_$V.log 2>&1; echo "$V dense exit $? $(date +%H:%M)"
$SQ -I driven/flux_metrics.py --run $R/flux_loc_$V --tag $V > results/metrics_$V.log 2>&1; echo "$V metrics exit $?"
$SQ -I driven/flux_xtalk_matrix.py --run $R/flux_loc_$V --tag $V > results/xtalk_$V.log 2>&1; echo "$V xtalk exit $?"
$SQ -I driven/flux_square_spectro.py --run $R/flux_loc_$V --line Q1 --amp 0.24 --tag $V --date $D > results/sq_Q1_$V.log 2>&1; echo "$V square exit $? $(date +%H:%M)"
echo B2 DONE
