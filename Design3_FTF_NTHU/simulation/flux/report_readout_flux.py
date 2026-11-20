"""Summarize Design3_FTF_NTHU electrical extractions against Design3_2 and Design3e."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I, glob_files as _glob_files, glob_paths as _glob_paths
_configure(__file__)
import json
from math import pi
from pathlib import Path
import sys

import pandas as pd
from scipy.constants import h, e

HERE = _P('Design3_FTF_NTHU')
RESULTS = _P('Design3_FTF_NTHU/results/general/legacy')
ROOT = _P('.')
OLD = _P('Design3_FTF_NTHU/results/general/legacy')
sys.path[:0] = [_I('Design3_FTF_NTHU'), _I('Design3_FTF_NTHU/simulation')]
from validate_gate import cap, extracted
from zz_analysis import fluxonium

NAMES = ('FX1', 'FX2')
PADS = {'FX1': ('pad_bot_FX1','pad_top_FX1'), 'FX2': ('pad_top_FX2','pad_bot_FX2')}
PHI0 = h/(2*e)


def read_json(path):
    return json.loads(path.read_text(encoding='utf8'))


def read_matrix(path):
    frame = pd.read_csv(path, index_col=0)
    return frame.loc[frame.index, frame.index]


def external_kappa(cap_ff, frequency, z=50.9):
    q = pi/(z*z*(cap_ff*1e-15*2*pi*frequency*1e9)**2)
    return frequency*1e3/q


def tee_cap(matrix, device):
    feed = cap.find_node(matrix, 'prime_cpw', 'feed')
    res = cap.find_node(matrix, f'second_cpw_tee_{device}', f'trace_res_{device}')
    value = -matrix.loc[res,feed]
    if value <= 0 or feed == res:
        raise ValueError(f'Invalid feedline/tee net for {device}')
    return float(value)


def mean_readout(device):
    slots = pd.read_csv(_relocate(RESULTS/f'readout_{device}_slots.csv'))
    fields = slots.assign(E2=slots.E**2).groupby(['mode','frequency_GHz','region']).E2.mean().unstack()
    fields['E_ratio'] = fields['near']/fields['far']
    choice = fields.sort_values('E_ratio').reset_index().iloc[-1]
    conv = pd.read_csv(_relocate(RESULTS/f'readout_{device}_convergence.csv'), index_col=0)
    mesh = pd.read_csv(_relocate(RESULTS/f'readout_{device}_mesh.csv'), index_col=0)
    i = int(choice['mode'])-1
    return dict(f_GHz=float(choice.frequency_GHz), mode=i+1, E_ratio=float(choice.E_ratio),
                last_shift_MHz=1000*float(conv.iloc[-1,i]-conv.iloc[-2,i]),
                last_max_delta_pct=float(mesh.iloc[-1,-1]), passes=len(conv))


def m_fixture(tag):
    summary = read_json(_relocate(RESULTS/f'm_{tag}_summary.json'))
    z = pd.read_csv(_relocate(RESULTS/f'm_{tag}_Z.csv'))
    conv = pd.read_csv(_relocate(RESULTS/f'm_{tag}_convergence.csv'), index_col=0)
    messages = (_relocate(RESULTS/f'm_{tag}_messages.log')).read_text(encoding='utf8')
    problem = [line for line in messages.splitlines() if ('[error]' in line.lower() or
               ('[warning]' in line.lower() and 'port' in line.lower()))]
    return {'M_pH': float(summary['M_pH_fmin']), 'I_phi0_mA': float(summary['I_phi0_mA']),
            'sweep_spread_pct': float(summary['M_spread_pct']),
            'reciprocity_error': float(z['Z21_Z12_mismatch'].max()),
            'last_delta_s': float(conv.iloc[-1,-1]), 'passes': len(conv), 'port_or_error_messages': problem}


def main():
    baseline = read_json(_P('Design3_FTF_NTHU/results/general/legacy/electrical_parameters.json'))
    old_tc = read_json(_P('Design3_FTF_NTHU/results/general/legacy/tc_68GHz_1pH_summary.json'))
    old_loaded = read_matrix(_P('Design3_FTF_NTHU/results/general/legacy/d32_matrix_cmat_fF.csv'))
    final_path = _P('Design3_FTF_NTHU/results/general/legacy/arms6p7_18p9_cmat_fF.csv')
    mesh_path = _P('Design3_FTF_NTHU/results/general/legacy/arms6p7_18p9_pass11_cmat_fF.csv')
    final = read_matrix(final_path)
    previous = read_matrix(mesh_path)
    ecj, _ = extracted(final_path)
    mesh_ecj, _ = extracted(mesh_path)
    feed = read_matrix(_P('Design3_FTF_NTHU/results/general/legacy/readout_kappa_cmat_fF.csv'))
    feed_passes = sorted(_glob_paths(RESULTS, 'readout_kappa_pass*_cmat_fF.csv'),
                         key=lambda path: int(path.stem.split('_')[2][4:]))
    outputs = dict(EcJ=ecj, pass11_EcJ=mesh_ecj, readout={}, mutual={})
    for device in NAMES:
        row = mean_readout(device)
        f = row['f_GHz']
        cf = tee_cap(feed,device)
        line = cap.find_node(final,f'readout_line_{device}',f'res_{device}')
        old_line = cap.find_node(old_loaded,f'readout_line_{device}',f'res_{device}')
        get_g = lambda matrix, node, frequency: abs(cap.readout_g(
            matrix,PADS[device],node,f_r=frequency*1e9,Z0=50.9,verbose=False)['g/2pi [MHz]'])
        _, n = fluxonium(ecj[f'Ec_{device} [GHz]'],4,.6)
        charge_element = float(abs(n[0,1]))
        old_f = baseline['eigen'][device]['f']
        old_g = get_g(old_loaded,old_line,old_f)
        _, old_n = fluxonium(baseline['EcJ'][f'Ec_{device} [GHz]'],4,.6)
        g = get_g(final,line,f)
        row.update(tee_C_fF=cf, kappa_MHz=external_kappa(cf,f),
                   kappa_pass_MHz={path.stem:external_kappa(tee_cap(read_matrix(path),device),f) for path in feed_passes},
                   qubit_res_per_n_MHz=g, g01_MHz=g*charge_element, n01=charge_element,
                   per_n_pass11_MHz=get_g(previous,line,f),
                   old_f_GHz=old_f, old_kappa_MHz=external_kappa(baseline['kappa'][device]['C_1_fF'],old_f),
                   old_per_n_MHz=old_g, old_g01_MHz=old_g*float(abs(old_n[0,1])))
        outputs['readout'][device] = row
    fixtures = {tag:m_fixture(tag) for tag in ('fx1','fx2','tc_left','tc_right')}
    outputs['mutual_fixtures'] = fixtures
    for q, tags in (('FX1',['fx1']),('FX2',['fx2']),('TC',['tc_left','tc_right'])):
        m = sum(abs(fixtures[tag]['M_pH']) for tag in tags)/len(tags)
        old_m = old_tc['mutual_mean_pH'] if q == 'TC' else baseline['mutual'][q]['M_pH']
        outputs['mutual'][q] = dict(M_pH=m, I_phi0_mA=PHI0*1e15/m,
                                    old_M_pH=old_m, old_I_phi0_mA=PHI0*1e15/old_m)
    (_P('Design3_FTF_NTHU/results/general/legacy/readout_flux_parameters.json')).write_text(json.dumps(outputs, indent=2), encoding='utf8')
    print(json.dumps(outputs, indent=2))


if __name__ == '__main__':
    main()
