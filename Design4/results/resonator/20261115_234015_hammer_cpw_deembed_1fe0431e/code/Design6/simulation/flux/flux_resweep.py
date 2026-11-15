"""Re-use a solved hfss_flux.py design (its adaptive mesh is in the AEDT project): add a frequency sweep on that mesh,
solve only the sweep and export the 6 x 6 Z matrix.

  --dense : Interpolating sweep 1 MHz - 3 GHz (for the pulse-distortion study, flux_pulse.py), saved as Z_dense.npy /
            f_dense.npy next to the run's Z.npy (which is not touched).
  default : short Discrete sweep (npts, 10 MHz - 1 GHz), saved as Z.npy / f.npy, then hfss_flux.analyse().

Usage:  python driven/flux_resweep.py --tag loc_fl10 --dense [--npts 600] [--design flux_loc_fl10_1556_hfss]
"""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np

HERE = _S('Design6/simulation/flux')
SIM = _S('Design6/simulation')
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', required=True)
    ap.add_argument('--design', default=None)
    ap.add_argument('--dense', action='store_true')
    ap.add_argument('--npts', type=int, default=None)
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, f'flux_{a.tag}'))
    os.makedirs(out, exist_ok=True)
    from pyEPR.ansys import HfssApp
    sys.path.insert(0, _I('Design5/simulation'))
    import save_aedt                     # activate the project the flux run used (env AEDT_PROJECT), not whatever is active
    save_aedt.use_project(os.environ.get('AEDT_PROJECT', _S('Design6/simulation/aedt/Design6.aedt')))
    proj = HfssApp().get_app_desktop().get_active_project()
    names = sorted(d.name for d in proj.get_designs() if d.name.startswith(f'flux_{a.tag}_'))
    dname = a.design or names[-1]
    des = [d for d in proj.get_designs() if d.name == dname][0]
    print('design', dname, flush=True)
    setup = des.get_setup('Setup')
    if a.dense:
        npts = a.npts or 600
        sw_name, kw = 'Dense', dict(start_ghz=0.001, stop_ghz=3.0, count=npts, type='Interpolating')
    else:
        npts = a.npts or 4
        sw_name, kw = f'Sweep{npts}', dict(start_ghz=0.01, stop_ghz=1.0, count=npts, type='Discrete')
    if sw_name not in setup.get_sweep_names():
        setup.insert_sweep(name=sw_name, save_fields=False, **kw)
    t0 = time.time()
    des._design.Analyze(f'Setup : {sw_name}')
    print(f'sweep {sw_name} solved in {(time.time() - t0) / 60:.1f} min', flush=True)
    sweep = setup.get_sweep(sw_name)
    fmts = [f'Z{i}{j}' for i in range(1, 7) for j in range(1, 7)]
    freqs, vals = sweep.get_network_data(fmts)
    freqs = np.asarray(freqs, float)
    Z = np.zeros((len(freqs), 6, 6), complex)
    for n, v in zip(fmts, vals):
        Z[:, int(n[1]) - 1, int(n[2]) - 1] = np.asarray(v)
    suffix = '_dense' if a.dense else ''
    np.save(_relocate(os.path.join(out, f'Z{suffix}.npy')), Z); np.save(_relocate(os.path.join(out, f'f{suffix}.npy')), freqs)
    json.dump(dict(design=dname, sweep=sw_name, sweep_options=kw, minutes=round((time.time() - t0) / 60, 1)),
              open(_relocate(os.path.join(out, f'run{suffix}.json')), 'w'))
    print('saved', out, len(freqs), 'points', flush=True)
    if not a.dense:
        import hfss_flux
        hfss_flux.analyse(out)


if __name__ == '__main__':
    main()
