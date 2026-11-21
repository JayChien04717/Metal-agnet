"""Compute fixed-LJ r9 qubit parameters and the bare-readout dispersive shift."""
from argparse import ArgumentParser
import json
from pathlib import Path
import sys
import traceback

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment

import cz_scqubits as C

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / 'code/qubit'), str(ROOT / 'code/general')]
import cap_matrix as cm
from run_record import create_run, finish_run

Q3D = ROOT / 'Design6_TXT/results/qubit/20261117_161936_q3d_core_r9_3b00f1f5'
CZ = ROOT / 'Design6_TXT/results/qubit/20261117_181834_cz_r9_20f32abc'


def dispersive(energy, charge, fr, gn):
    n = len(energy)
    a = np.diag(np.sqrt(np.arange(1, n)), 1)
    H = (np.kron(np.diag(energy), np.eye(n))
         + np.kron(np.eye(n), np.diag(np.arange(n) * fr))
         + gn * np.kron(charge, a + a.T))
    eig, vectors = eigh(H)
    states = [0, 1, n, n + 1]  # |0,0>, |0,1>, |1,0>, |1,1>
    overlap = abs(vectors[states]) ** 2
    rows, cols = linear_sum_assignment(-overlap)
    e00, e01, e10, e11 = eig[cols]
    return dict(two_chi_MHz=float((e11 - e10 - e01 + e00) * 1000),
                minimum_overlap=float(np.min(overlap[rows, cols])))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('s21_comparison', type=Path)
    args = parser.parse_args()
    s21 = args.s21_comparison.resolve()
    params = CZ / 'params_r9_q3d.json'
    cmat = Q3D / 'raw/cmat_fF.csv'
    folder = create_run('Design6_TXT', 'qubit', 'readout_qubit_r9', entrypoint=__file__,
                        inputs=[params, cmat, Q3D / 'summary_out.json', s21 / 'summary.json'],
                        upstream=[str(Q3D), str(CZ), str(s21)],
                        parameters=dict(EM_rerun=False, geometry_revision='2026-11-16_wafer_r9',
                                        EJ='layout LJ fixed', charge_cutoffs=[15, 25],
                                        local_levels_and_photon_cutoffs=[6, 8], XY_load_ohm=50,
                                        chi_model='isolated Xmon plus its bare resonator, full charge interaction'))
    print('RUN_FOLDER=' + str(folder), flush=True)
    try:
        p = C.load_params(params)
        core = json.loads((Q3D / 'summary_out.json').read_text(encoding='utf-8'))
        loaded = json.loads((s21 / 'summary.json').read_text(encoding='utf-8'))['channels']
        capacitances = pd.read_csv(cmat, index_col=0)
        result = {}
        for q in ('Q1', 'Q2', 'TC'):
            fits = []
            for n, ncut in ((6, 15), (8, 25)):
                qubit = C.scq.Transmon(EJ=p['EJsum_GHz'][q], EC=p['Ec_GHz'][q], ng=0, ncut=ncut)
                vals, vectors = qubit.eigensys(evals_count=n)
                energy = vals - vals[0]
                charge = vectors.conj().T @ qubit.n_operator() @ vectors
                n01 = float(abs(charge[0, 1]))
                record = dict(f01_GHz=float(energy[1]), alpha_MHz=float((energy[2] - 2 * energy[1]) * 1000),
                              n01=n01, local_levels=n, charge_cutoff=ncut)
                if q != 'TC':
                    fr = loaded[q]['candidate']['bare_fr_GHz']['hfss']
                    cap = core['final']['readout'][q]
                    gn = cap['g/2pi [MHz]'] * fr / cap['f_r [GHz]'] / 1000
                    record.update(fr_GHz=fr, g_charge_MHz=gn * 1000, g01_MHz=gn * n01 * 1000,
                                  **dispersive(energy, charge, fr, gn))
                    line = 'XY1' if q == 'Q1' else 'XY2'
                    nmap = core['net_mapping']
                    record['xy'] = cm.charge_line(
                        capacitances, (nmap[q],), nmap[line], energy[1] * 1e9, n01)
                fits.append(record)
            result[q] = dict(main=fits[1], truncation_check=fits[0],
                             Lj_nH={k: v for k, v in p['Lj_nH'].items() if k.startswith(q)},
                             EJsum_GHz=p['EJsum_GHz'][q], Ec_GHz=p['Ec_GHz'][q])
        finish_run(folder, status='completed', summary=result)
        print(json.dumps(result, indent=2, default=float))
    except Exception:
        (folder / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        finish_run(folder, status='failed', summary={'error': traceback.format_exc()})
        raise


if __name__ == '__main__':
    main()
