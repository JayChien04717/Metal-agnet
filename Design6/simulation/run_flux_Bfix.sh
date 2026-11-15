#!/bin/bash
# Re-run the dense sweep + metrics of variant B (its first dense sweep looked in the wrong AEDT project), after D,
# then the A-D comparison. Arg: output file of the D run to wait for.
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
export AEDT_PROJECT="$(pwd -W)/aedt/Design6_flux.aedt"
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
SQ=/c/Users/zx123/miniforge3/envs/scqenv/python.exe
R=results/2026-11-14
until grep -q "D DONE" "$1"; do sleep 30; done
V=Lshort
$PY -I driven/flux_resweep.py --tag loc_$V --dense > results/flux_dense_$V.log 2>&1; echo "$V dense exit $? $(date +%H:%M)"
$SQ -I driven/flux_metrics.py --run $R/flux_loc_$V --tag $V > results/metrics_$V.log 2>&1; echo "$V metrics exit $?"
$SQ -I driven/flux_xtalk_matrix.py --run $R/flux_loc_$V --tag $V > results/xtalk_$V.log 2>&1; echo "$V xtalk exit $?"
$SQ -I driven/flux_square_spectro.py --run $R/flux_loc_$V --line Q1 --amp 0.24 --tag $V > results/sq_Q1_$V.log 2>&1; echo "$V square exit $?"
$PY -I driven/compare_flux_variants.py > results/compare_flux.log 2>&1; echo "compare exit $?"
echo BFIX DONE
