# Wait until AEDT has finished the eigenmode solve of v6_full_eig_hfss (no hf3d for 3 min), then the driven design
# in the same project (Project1) with the eigenmode mesh imported, ONE pass on that mesh (no refinement), Fast sweep,
# fit. Log: results\v6_full.log
$py = 'C:\Users\QEL\miniforge-pypy3\envs\qmetal\python.exe'
Set-Location $PSScriptRoot
$idle = 0
while ($idle -lt 6) {
    if (Get-Process hf3d -ErrorAction SilentlyContinue) { $idle = 0 } else { $idle++ }
    Start-Sleep -Seconds 30
}
"eigenmode solve finished $(Get-Date -Format HH:mm)"
$env:SKIP_EIG = '1'; $env:CONTINUE_EIG = ''; $env:EIG_DESIGN = 'v6_full_eig_hfss'; $env:DRV_MAXP = '1'
$env:AEDT_PROJECT = (Join-Path (Split-Path $PSScriptRoot) 'Project1.aedt')
$env:EIG_NOTE = 'run_all -MeshLink (1st order, 10/5 um initial mesh); Delta F changed 0.05 -> 0.1 % in AEDT by Jay after pass 3; continued in AEDT to the end'
$env:RUN_NOTE = 'Design6 verification, full channel (no C_in), eigenmode mesh (1st order, Delta F 0.1 %) imported into driven 2nd order, 1 pass (no refinement), short feed, HFSS wire bonds'
cmd /c "`"$py`" -W ignore driven\hfss6_meshlink.py v6_full > results\v6_full.log 2>&1"
"driven (linked) exit $LASTEXITCODE"
