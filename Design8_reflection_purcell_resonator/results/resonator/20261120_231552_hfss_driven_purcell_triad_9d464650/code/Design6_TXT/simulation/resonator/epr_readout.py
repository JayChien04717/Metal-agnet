"""pyEPR on an already solved Design6_TXT readout-channel eigenmode design (hfss_readout.py), no re-solve.

The qubit SQUID = two lumped junctions JJ_rect_Lj_<Q>_rect_jj1 / _jj2 (variables Lj_<Q>a / Lj_<Q>b).
Output: <run folder>/epr.json (f_ND, chi_ND = numerical-diagonalisation cross-Kerr matrix, MHz), then summary.json.

Usage: python epr_readout.py Q1 [--tag v8] [--design roQ1_v8_1115_hfss]   (default: newest roQ<n>_<tag>_* design)
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

HERE = _S('Design6_TXT/simulation/resonator')
SIM = _S('Design6_TXT/simulation')
sys.path[:0] = [HERE]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('qubit')
    ap.add_argument('--tag', default='v8')
    ap.add_argument('--design', default=None)
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    ap.add_argument('--npz', default=None, help='pyEPR data file of a finished do_EPR_analysis (offline, no AEDT)')
    a = ap.parse_args()
    q = a.qubit
    out = _relocate(os.path.join(SIM, 'results', a.date, f'hfss_readout_{q}_{a.tag}'))
    from pyEPR import ProjectInfo, DistributedAnalysis, QuantumAnalysis
    if a.npz:
        data, dname = a.npz, os.path.basename(os.path.dirname(a.npz))
    else:
        from pyEPR.ansys import HfssApp
        proj = HfssApp().get_app_desktop().get_active_project()
        names = sorted(d.name for d in proj.get_designs() if d.name.startswith(f'ro{q}_{a.tag}_'))
        dname = a.design or names[-1]
        print('design', dname, flush=True)
        pinfo = ProjectInfo(project_path=_S('Design6_TXT/simulation/aedt'), project_name='Design6_TXT', design_name=dname)
        for k, ab in ((1, 'a'), (2, 'b')):
            pinfo.junctions[f'j{q}{ab}'] = dict(Lj_variable=f'Lj_{q}{ab}', rect=f'JJ_rect_Lj_{q}_rect_jj{k}',
                                                line=f'JJ_Lj_{q}_rect_jj{k}_', Cj_variable='Cj')
        pinfo.validate_junction_info()
        epra = DistributedAnalysis(pinfo)
        epra.do_EPR_analysis()
        data = epra.data_filename
    qa = QuantumAnalysis(data)
    r = qa.analyze_all_variations(cos_trunc=8, fock_trunc=7, print_result=False)
    rr = r[list(r.keys())[0]]
    f_nd = np.array(rr['f_ND'], float)
    f_nd = f_nd / 1e3 if f_nd.max() > 100 else f_nd
    f_0 = np.array(rr['f_0'], float)
    f_0 = f_0 / 1e3 if f_0.max() > 100 else f_0
    chi = np.array(rr['chi_ND'], float)
    epr = dict(design=dname, epr_data=data, f_ND_GHz=f_nd.tolist(), f_0_GHz=f_0.tolist(), chi_ND_MHz=chi.tolist(),
               note='pyEPR chi_ND: diagonal = anharmonicity alpha (MHz), off-diagonal = full dispersive shift 2chi')
    json.dump(epr, open(_relocate(os.path.join(out, 'epr.json')), 'w'), indent=1)
    print('f_ND (GHz):', np.round(f_nd, 4), '\nchi_ND (MHz):\n', np.round(chi, 4), flush=True)
    import hfss_readout
    hfss_readout.summarize(out)


if __name__ == '__main__':
    main()
