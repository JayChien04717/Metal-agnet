# Design6 verification chain: Q3D (full, filter) -> [HFSS filter] -> HFSS full. Logs: results/<tag>.log (moved with the run).
# Usage: run_all.ps1 [-SkipQ3D] [-WithFilter]
#   The filter-only HFSS run is optional: with no C_in the full-channel S21 alone fixes dL_filter, dL_res, C_p, C_j.
$ErrorActionPreference = 'Continue'
$py = 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe'
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force results | Out-Null
# one AEDT session, one project per Design (simulation\aedt\Design6.aedt, opened by the scripts): nothing is closed
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
if ($args -contains '-MeshLink') {
    # mesh adapted by an eigenmode solve, imported into the driven design
    $env:RUN_NOTE = "Design6 verification, full channel (no C_in), eigenmode mesh -> driven (mesh link), short feed, HFSS wire bonds, initial mesh 10 um lines / 5 um capacitors"
    cmd /c "`"$py`" -W ignore driven\hfss6_meshlink.py v6_full > results\v6_full.log 2>&1"
    "hfss full (mesh link) exit $LASTEXITCODE"
} else {
    $env:RUN_NOTE = "Design6 verification, full channel (no C_in), $mesh"
    cmd /c "`"$py`" -W ignore driven\hfss6.py full v6_full > results\v6_full.log 2>&1"
    "hfss full exit $LASTEXITCODE"
}
