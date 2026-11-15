#!/bin/bash
# Design6 HFSS queue (one AEDT, one job at a time): readout channels, flux-line variants. Started after hfss_core.
PY=/c/Users/zx123/miniforge3/envs/qmetal/python.exe
export PYTHONIOENCODING=utf-8
cd "$(dirname "$0")"
# (core sweep already done)\|Traceback" results/hfss_core_v8.log; do sleep 20; done
echo "core done $(date +%H:%M)"
$PY -I eigenmode/hfss_readout.py Q1 --tag v8s > results/hfss_readout_Q1_v8s.log 2>&1; echo "readout Q1 exit $? $(date +%H:%M)"
$PY -I eigenmode/hfss_readout.py Q2 --tag v8s > results/hfss_readout_Q2_v8s.log 2>&1; echo "readout Q2 exit $? $(date +%H:%M)"
OV0='{"coupler": {"squid_contact_overlap": "0um", "squid_shift_toward_stem": "0um"}}'
$PY -I driven/hfss_flux.py --tag v8_fl10 --overrides "$OV0" > results/hfss_flux_v8_fl10.log 2>&1; echo "flux V0 exit $? $(date +%H:%M)"
OV1='{"coupler": {"squid_contact_overlap": "0um", "squid_shift_toward_stem": "0um", "flux_bias_line_options": {"make_fbl": true, "fbl_sep": "5um", "fbl_length": "60um", "fbl_x": "-25um"}}, "qflux": {"fbl_sep": "5um", "fbl_length": "60um"}}'
$PY -I driven/hfss_flux.py --tag v8_fl5 --overrides "$OV1" > results/hfss_flux_v8_fl5.log 2>&1; echo "flux V1 exit $? $(date +%H:%M)"
