# Design6 verification chain: Q3D (full, filter) -> [HFSS filter] -> HFSS full. Logs: results/<tag>.log (moved with the run).
# Usage: run_all.ps1 [-SkipQ3D] [-WithFilter]
#   The filter-only HFSS run is optional: with no C_in the full-channel S21 alone fixes dL_filter, dL_res, C_p, C_j.
$ErrorActionPreference = 'Continue'
$py = 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe'
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force results | Out-Null
# close projects left open in AEDT (an earlier run, already saved or stopped)
& $py -W ignore -c "from pyEPR.ansys import HfssApp; d = HfssApp().get_app_desktop(); [d._desktop.CloseProject(p.name) for p in d.get_projects()]; print('AEDT projects closed')"
if ($args -notcontains '-SkipQ3D') { foreach ($v in 'full', 'filter') {
    cmd /c "`"$py`" -W ignore q3d\q3d6.py $v > results\q3d6_$v.log 2>&1"
    "q3d6 $v exit $LASTEXITCODE"
} }
$mesh = 'short feed, HFSS wire bonds, initial mesh 10 um lines / 5 um capacitors, 3-15 passes, 20 % refinement'
if ($args -contains '-WithFilter') {
    $env:RUN_NOTE = "Design6 verification, filter only (C_p + filter, no C_in), $mesh"
    cmd /c "`"$py`" -W ignore driven\hfss6.py filter v6_filter > results\v6_filter.log 2>&1"
    "hfss filter exit $LASTEXITCODE"
}
$env:RUN_NOTE = "Design6 verification, full channel (no C_in), $mesh"
cmd /c "`"$py`" -W ignore driven\hfss6.py full v6_full > results\v6_full.log 2>&1"
"hfss full exit $LASTEXITCODE"
