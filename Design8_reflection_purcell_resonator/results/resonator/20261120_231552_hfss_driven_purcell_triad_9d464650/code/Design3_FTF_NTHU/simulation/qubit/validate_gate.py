"""(Copied from Design3e so Design3_FTF_NTHU runs on its own.) Verify a Q3D candidate with the previous Design3_2 multilevel CZ model."""

# Project paths: explicit migration map; solver/geometry parameters are unchanged.
from pathlib import Path as _ProjectPath
import sys as _project_sys
_project_root = _ProjectPath(__file__).resolve().parents[3]
_project_sys.path.insert(0, str(_project_root / 'code' / 'general'))
from project_paths import path as _P, string_path as _S, relocate as _relocate, configure as _configure, import_path as _I
_configure(__file__)
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = _P('Design3_FTF_NTHU')
ROOT = _P('.')
sys.path.insert(0, _I('note'))
import cap_matrix as cap

def extracted(path, allow_local=False):
    c = pd.read_csv(path, index_col=0)
    c = c.loc[c.index, c.index]
    raw = c.to_numpy()
    if not np.allclose(raw, raw.T, rtol=1e-5, atol=1e-6) or np.linalg.eigvalsh(raw).min() <= 0:
        raise ValueError('Q3D matrix must be symmetric and positive definite')
    devices = ('FX1', 'TC', 'FX2')
    readouts = [cap.find_node(c, f'readout_line_{q}', f'res_{q}') for q in ('FX1','FX2')]
    readouts += [q for q in ('readout_line_TC','res_TC') if q in c.index][:1]
    nodes = cap.DEFAULT_NODES[:5] + readouts
    if not allow_local and min(c.loc[r,r] for r in readouts) < 500:
        raise ValueError('Use the full-resonator Q3D result')
    electrical = cap.analyse(c, nodes, verbose=False)
    lines = [f'charge_line_{q}' for q in devices]
    matrix = cap.reduced_matrix(c, nodes+lines)
    transform = np.linalg.inv(cap.transform_matrix(len(matrix)))
    inverse = np.linalg.inv(transform.T @ matrix @ transform)
    indices = list(range(len(nodes),len(nodes)+len(lines)))
    beta = inverse[np.ix_([1,2,4],indices)] / np.diag(inverse)[indices]
    return electrical, beta

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix', type=Path)
    parser.add_argument('--model', type=Path, default=_P('Design3_FTF_NTHU/simulation/qubit'))
    parser.add_argument('--allow-local',action='store_true')
    parser.add_argument('--ramp-points', type=int, default=240)
    parser.add_argument('--fx-levels', type=int, default=8)
    parser.add_argument('--tc-levels', type=int, default=6)
    parser.add_argument('--gate-ns', type=int, nargs='+', default=[25, 31, 35, 45, 55, 65, 80, 95])
    parser.add_argument('--ramp-ns', type=int, nargs='+', default=[2,3,4,5,7])
    parser.add_argument('--basis-candidates', type=int, default=8)
    args = parser.parse_args()
    electrical, beta = extracted(args.matrix,args.allow_local)
    sys.path.insert(0, str(args.model.resolve()))
    import verify as ramp
    import sweep as model
    from sequence import ramp_U, flat_U
    on, off = ramp.TC_CZ_GHZ, ramp.TC_IDLE_GHZ
    model.ec1,model.ect,model.ec2,model.j1,model.j2,model.j12 = electrical.values()
    model.beta = beta
    nf, nt = args.fx_levels, args.tc_levels
    w, _, ids, ops = model.setup(on, nf, nt)
    path = ramp.ramp_path(off, on, args.ramp_points, nf, nt)
    original_setup = model.setup
    model.setup = lambda f: original_setup(f, nf, nt)
    try:
        cases = model.candidates(on)
    finally:
        model.setup = original_setup
    output = _relocate(_P('Design3_FTF_NTHU/results/general/legacy')/f'{args.matrix.stem}_gate_{args.ramp_points}steps')
    output.with_suffix('.json').write_text(json.dumps(dict(electrical=electrical,beta=beta.tolist(),
        idle=model.metrics(off,nf,nt),on=model.metrics(on,nf,nt),EJ_GHz=path[-2:]),indent=2))
    candidates = sorted((row for row, *_ in cases if row['line']=='TC'),
                        key=lambda row: row['delta_sel_MHz'], reverse=True)[:args.basis_candidates]
    X = np.einsum('q,qij->ij', beta[:,1], ops)
    idle_ids = path[2]
    rows = []
    for tr in args.ramp_ns:
        before = ramp_U(path,tr)[:,idle_ids]
        after = ramp_U(path,tr,True)[idle_ids,:]
        for gate in args.gate_ns:
            T = gate-2*tr
            if T < 14:
                continue
            for case in candidates:
                i,j = case['i'],case['j']
                if abs(X[j,i]) < 1e-7:
                    continue
                carrier = w[j]-w[i]
                amp0 = 2/(T*abs(X[j,i]))
                def score(x):
                    U = after @ flat_U(w,X,carrier+x[1]/T,T,amp0*x[0]) @ before
                    phases = np.angle(np.diag(U))
                    corrected = np.exp(-1j*np.array([phases[0],phases[1],phases[2],phases[1]+phases[2]-phases[0]]))
                    trace = np.dot(corrected*np.array([1,1,1,-1]),np.diag(U))
                    return 1-(np.sum(abs(U)**2)+abs(trace)**2)/20
                result = minimize(score,[1,0],method='Nelder-Mead',bounds=[(.65,1.4),(-.55,.55)],
                                  options=dict(initial_simplex=[[1,0],[1.04,0],[1,.05]],maxiter=65))
                rows.append(dict(**case,total_ns=gate,mw_ns=T,ramp_ns=tr,error=result.fun,
                                 drive_GHz=carrier+result.x[1]/T,drive_scale_GHz=amp0*result.x[0],
                                 drive_peak_mV=amp0*result.x[0]*4.135667696e-3/2))
    with output.with_suffix('.csv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('Best at 31 ns', min((row for row in rows if row['total_ns']==31),key=lambda row:row['error'],default=None))
    print('Fastest F >= 99.9%', min((row for row in rows if row['error']<.001),key=lambda row:row['total_ns'],default=None))

if __name__=='__main__':
    main()
