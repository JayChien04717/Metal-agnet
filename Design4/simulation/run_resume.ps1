# Resume v6_full: continue the eigenmode design already open in AEDT (v6_full_eig_hfss) to its Delta F (0.1 %, set in
# AEDT by Jay), then the driven design with the imported mesh, the sweep and the fit. Log: results\v6_full.log
$py = 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe'
Set-Location $PSScriptRoot
$env:SKIP_EIG = '1'; $env:CONTINUE_EIG = '1'; $env:EIG_DESIGN = 'v6_full_eig_hfss'
$env:EIG_NOTE = 'started by run_all -MeshLink (1st order, 10/5 um initial mesh), stopped after pass 3; Delta F changed 0.05 -> 0.1 % in AEDT by Jay; continued from pass 3'
$env:RUN_NOTE = 'Design6 verification, full channel (no C_in), eigenmode mesh (1st order, Delta F 0.1 %) -> driven 2nd order (mesh link), short feed, HFSS wire bonds'
cmd /c "`"$py`" -W ignore driven\hfss6_meshlink.py v6_full > results\v6_full.log 2>&1"
"resume exit $LASTEXITCODE"
