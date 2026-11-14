from pathlib import Path
import os
import sys

design_dir = Path(__file__).resolve().parents[1]
root = design_dir.parent
out = design_dir / 'simulation/results/layout3_full'
out.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(design_dir / 'pattern design'))

import json
import numpy as np,pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import least_squares
sys.path.insert(0,str(root/'note'))
import cap_matrix as cap
matrix=pd.read_csv(out/'capacitance_fF.csv',index_col=0)
def node(name):
    return cap.find_node(matrix,'metal_'+name)
nodes=[node(n) for n in ('fx_left_bottom','fx_left_top','transmon',
                       'fx_right_top','fx_right_bottom','res_left','res_right')]
result=cap.analyse(matrix,nodes,verbose=False)
result['C_sigma_fF']=cap.total_capacitances(matrix,nodes)
result['charge_direct_fF']={line:{name:float(-matrix.loc[node(line),node(name)]) for name in
    ('fx_left_top','fx_left_bottom','transmon','fx_right_top','fx_right_bottom')}
    for line in ('charge_left','charge_right')}
conv=pd.read_csv(out/'q3d_convergence.csv')
result['q3d_passes']=len(conv)
result['q3d_final_delta_pct']=float(conv['Delta %'].iloc[-1])
for tag in ('readout', 'linear_readout', 'linear_readout_window1', 'linear_readout_window2'):
    sfile = out / (f'{tag}.csv' if 'window' in tag else f'{tag}_sparams.csv')
    if not sfile.exists():
        continue
    s = pd.read_csv(sfile)
    f = s.frequency_Hz.to_numpy() / 1e9
    z = np.column_stack((s.S21_re + 1j*s.S21_im, s.S11_re + 1j*s.S11_im))
    db = 20*np.log10(np.maximum(abs(z), 1e-15))
    minima, _ = find_peaks(-db[:, 0], prominence=.02)
    spacing_MHz = (f[1]-f[0])*1000
    fits = []
    for i in minima:
        center = f[i]
        neighbors = [abs(f[k]-center) for k in minima if k != i]
        half = min(.06, min(neighbors)*.45) if neighbors else .06
        select = abs(f-center) < half
        x = (f[select]-center)*1000
        y = z[select]
        if len(y) < 6:
            continue

        def design_matrix(df, width):
            return np.column_stack((np.ones_like(x), x, 1/(1 + 2j*(x-df)/width)))

        def err(p):
            design = design_matrix(*p)
            coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
            delta = design @ coefficients - y
            return np.r_[delta.real.ravel(), delta.imag.ravel()]

        shift = max(spacing_MHz, 1)
        opt = least_squares(err, [0, max(1, spacing_MHz)],
                            bounds=([-shift,.02],[shift,100]))
        df,width = opt.x
        fits.append(dict(f_GHz=center+df/1000, loaded_kappa_MHz=width,
            raw_dip_dB=db[i,0], samples_per_fwhm=width/spacing_MHz,
            fit_relative_residual=float(np.linalg.norm(opt.fun)/np.linalg.norm(y)),
            mode_assignment='unassigned', candidate_only=True))
    result[f'{tag}_candidates'] = fits
    result[f'{tag}_max_power_balance_error'] = float(np.max(abs(np.sum(abs(z)**2, axis=1)-1)))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9,4.5), layout='constrained')
    ax.plot(f, db[:,0], label='S21', lw=1)
    ax.plot(f, db[:,1], label='S11', lw=1)
    for entry in fits:
        ax.axvline(entry['f_GHz'], color='#666', lw=.6, ls=':')
    ax.set(xlabel='Frequency (GHz)', ylabel='Magnitude (dB)')
    ax.legend()
    fig.savefig(out/f'{tag}.png', dpi=160)
    plt.close(fig)

freqs=[5.9,6.1]

result['readout_frequency_basis'] = 'nominal 5.9 and 6.1 GHz, pending mode assignment'
result['readout']={}
result['all_readout_g_MHz_per_n']={}
for side,freq in zip(('left','right'),freqs):
    rn=node('res_'+side)
    result['readout'][side]=cap.readout_g(matrix,
        (node('fx_'+side+'_top'),node('fx_'+side+'_bottom')),rn,Z0=45.1,f_r=freq*1e9,verbose=False)
    pairs = {'FX1': (node('fx_left_top'), node('fx_left_bottom')),
             'TC': (node('transmon'),),
             'FX2': (node('fx_right_top'), node('fx_right_bottom'))}
    result['all_readout_g_MHz_per_n'][side] = {
        name: cap.readout_g(matrix, pads, rn, Z0=45.1, f_r=freq*1e9, verbose=False)['g/2pi [MHz]']
        for name, pads in pairs.items()
    }
    c=-matrix.loc[node('feed'),rn]
    result['readout'][side]['C_feed_fF']=float(c)
    result['readout'][side]['kappa_cap_model_MHz']=float(freq*1000/np.pi*(45.1*2*np.pi*freq*1e9*c*1e-15)**2)
(out/'electrical.json').write_text(json.dumps(result,indent=2),encoding='utf8')
print(json.dumps(result,indent=2))
