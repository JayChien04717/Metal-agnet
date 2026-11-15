#!/bin/bash
# Flux-line variant D (U hairpin: sep 9 um, legs 40 um apart, return leg shorted 120 um below the bar), same pipeline as
# run_flux_BC.sh. Waits until the B / C run has finished (single AEDT session).
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
export AEDT_PROJECT="$(pwd -W)/aedt/Design6_flux.aedt"
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
SQ=/c/Users/zx123/miniforge3/envs/scqenv/python.exe
R=results/2026-11-14
WAIT=${1:-}
if [ -n "$WAIT" ]; then until grep -q "ALL DONE" "$WAIT"; do sleep 30; done; fi
V=U
OV='{"coupler": {"flux_bias_line_options": {"make_fbl": true, "fbl_sep": "9um", "fbl_length": "40um", "fbl_x": "0um", "fbl_type": "U", "fbl_drop": "120um"}}, "qflux": {"fbl_sep": "9um", "fbl_length": "40um", "fbl_x": "0um", "fbl_type": "U", "fbl_drop": "120um"}}'
$PY -I driven/hfss_flux.py --tag loc_$V --overrides "$OV" > results/hfss_flux_loc_$V.log 2>&1; echo "$V hfss exit $? $(date +%H:%M)"
$PY -I driven/flux_resweep.py --tag loc_$V --dense > results/flux_dense_$V.log 2>&1; echo "$V dense exit $? $(date +%H:%M)"
$SQ -I driven/flux_metrics.py --run $R/flux_loc_$V --tag $V > results/metrics_$V.log 2>&1; echo "$V metrics exit $?"
$SQ -I driven/flux_xtalk_matrix.py --run $R/flux_loc_$V --tag $V > results/xtalk_$V.log 2>&1; echo "$V xtalk exit $?"
$SQ -I driven/flux_square_spectro.py --run $R/flux_loc_$V --line Q1 --amp 0.24 --tag $V > results/sq_Q1_$V.log 2>&1; echo "$V square exit $? $(date +%H:%M)"
echo D DONE
