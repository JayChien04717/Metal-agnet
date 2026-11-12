"""Add an Interpolating sweep to an already solved driven design (same mesh) and save its S21.

The Fast sweep can lose a narrow, weakly driven resonance; the Interpolating sweep solves the full field at points it
places itself (adds points where S changes fast), so a resonator line is not lost.

Usage: python driven/interp_sweep.py <design name> <out folder> [--f0 6.47 --f1 6.50 --step 0.02 --maxsol 60]
Output: <out>/S_interp.csv
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('design'); ap.add_argument('out')
ap.add_argument('--f0', type=float, default=6.47); ap.add_argument('--f1', type=float, default=6.50)
ap.add_argument('--step', type=float, default=0.02, help='output step (MHz)')
ap.add_argument('--maxsol', type=int, default=60, help='max field solutions of the interpolating sweep')
ap.add_argument('--tol', type=float, default=0.2, help='interpolation tolerance (%%)')
a = ap.parse_args(); a.out = os.path.abspath(a.out)

SIM = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(SIM)), 'Design5', 'simulation'))
import save_aedt
save_aedt.use_project(save_aedt.project_file(SIM))
from pyEPR.ansys import HfssApp
proj = HfssApp().get_app_desktop().get_active_project()
des = proj._project.SetActiveDesign(a.design)
mod = des.GetModule('AnalysisSetup')
name = f'Interp_{time.strftime("%H%M")}'
count = int(round((a.f1 - a.f0) * 1e3 / a.step)) + 1
mod.InsertFrequencySweep('Setup', [
    f'NAME:{name}', 'IsEnabled:=', True, 'RangeType:=', 'LinearCount', 'RangeStart:=', f'{a.f0}GHz',
    'RangeEnd:=', f'{a.f1}GHz', 'RangeCount:=', count, 'Type:=', 'Interpolating', 'SaveFields:=', False,
    'SaveRadFields:=', False, 'InterpTolerance:=', a.tol, 'InterpMaxSolns:=', a.maxsol, 'InterpMinSolns:=', 0,
    'InterpMinSubranges:=', 1, 'ExtrapToDC:=', False])
t0 = time.time()
print(f'{a.design}: interpolating sweep {name} {a.f0}-{a.f1} GHz ({count} pts, tol {a.tol} %, <= {a.maxsol} solutions)', flush=True)
des.Analyze(f'Setup : {name}')
print(f'solved in {(time.time() - t0) / 60:.1f} min', flush=True)
rep = des.GetModule('ReportSetup')
csv = os.path.join(a.out, 'S_interp_raw.csv').replace(os.sep, '/')
rname = f'S21_{name}'
rep.CreateReport(rname, 'Modal Solution Data', 'Data Table', f'Setup : {name}', ['Domain:=', 'Sweep'],
                 ['Freq:=', ['All']], ['X Component:=', 'Freq', 'Y Component:=', ['re(S(1,2))', 'im(S(1,2))']])
rep.ExportToFile(rname, csv)
df = pd.read_csv(csv)
f = df.iloc[:, 0].values.astype(float)
f = f / 1e9 if f.max() > 1e3 else (f / 1e3 if f.max() > 100 else f)       # Hz / MHz -> GHz
out = pd.DataFrame({'S21': [f'({r}{i:+}j)' for r, i in zip(df.iloc[:, 1], df.iloc[:, 2])]}, index=pd.Index(f, name='f_GHz'))
out.to_csv(os.path.join(a.out, 'S_interp.csv'))
proj._project.Save()
print('saved', os.path.join(a.out, 'S_interp.csv'), len(f), 'points', flush=True)
