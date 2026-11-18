"""Design3 flux-line mutual inductance: converged 10-100 MHz run + broadband run, then the figure."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import os, sys as _sys; _sys.path[:0] = [_I('Design3/simulation')]
from qiskit_metal import Dict
import hfss_mutual_fx as m
import plot_mutual_z

m.run('mutual_d3')
m.run('mutual_d3_broadband', setup=dict(sweep_setup=Dict(name='Sweep', start_ghz=0.01, stop_ghz=6.0, count=300,
                                                         step_ghz=None, type='Interpolating', save_fields=False)))
plot_mutual_z.main('mutual_d3_broadband', 'mutual_d3')
