# Design6 verification chain: Q3D (full, filter) -> HFSS filter -> HFSS full. Logs: results/<tag>.log (moved with the run).
$ErrorActionPreference = 'Continue'
$py = 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe'
Set-Location $PSScriptRoot
New-Item -ItemType Directory -Force results | Out-Null
# close the project left open in AEDT (chip6_v1, already saved)
& $py -W ignore -c "from pyEPR.ansys import HfssApp; d = HfssApp().get_app_desktop(); [d._desktop.CloseProject(p.name) for p in d.get_projects()]; print('AEDT projects closed')"
if ($args -notcontains '-SkipQ3D') { foreach ($v in 'full', 'filter') {
    cmd /c "`"$py`" -W ignore q3d6.py $v > results\q3d6_$v.log 2>&1"
    "q3d6 $v exit $LASTEXITCODE"
} }
$env:RUN_NOTE = 'Design6 verification, filter only (C_p + filter, no C_in), short feed, HFSS wire bonds, initial mesh 5 um on CPW traces'
cmd /c "`"$py`" -W ignore hfss6.py filter v6_filter > results\v6_filter.log 2>&1"
"hfss filter exit $LASTEXITCODE"
$env:RUN_NOTE = 'Design6 verification, full channel (no C_in), short feed, HFSS wire bonds, initial mesh 5 um on CPW traces'
cmd /c "`"$py`" -W ignore hfss6.py full v6_full > results\v6_full.log 2>&1"
"hfss full exit $LASTEXITCODE"
