"""Narrow Discrete S21 sweeps around the eigenmodes of the full chain (eig_full) to look for
the narrow hanger-resonator dips (model: 0.2-0.4 MHz wide, much narrower than earlier 2 MHz grids).

Usage:  python hfss_dips.py [f1 f2 ...]   (GHz; default = the four modes of eig_full in 6.5-6.95 GHz)
"""
import os
import sys
import warnings

warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design'))
import hfss_combined as hc  # noqa: E402

HALF_MHZ, STEP_MHZ = 10.0, 0.25


def main(freqs):
    rows = []
    for f0 in freqs:
        hc.SWEEP.update(start_ghz=f0 - HALF_MHZ * 1e-3, stop_ghz=f0 + HALF_MHZ * 1e-3,
                        count=int(round(2 * HALF_MHZ / STEP_MHZ)) + 1, type='Discrete')
        hc.SETUP.update(freq_ghz=f0, max_delta_s=0.01, max_passes=10)  # adapt mesh AT the resonance
        tag = f'dip3_{f0:.4f}'
        hc.run(tag)
        S = pd.read_csv(os.path.join(hc.RESULTS, f'{tag}_S.csv'), index_col=0)
        f = S.index.values.astype(float)
        s21 = 20 * np.log10(np.abs(np.array([complex(str(v).replace(' ', '')) for v in S['S21']])))
        i = int(np.argmin(s21))
        rows.append(dict(window_GHz=f0, f_min_GHz=f[i], S21_min_dB=s21[i], S21_max_dB=s21.max(),
                         depth_dB=s21.max() - s21[i]))
        print(f'WINDOW {f0:.4f}: min S21 {s21[i]:.2f} dB at {f[i]:.4f} GHz, window max {s21.max():.2f} dB')
    pd.DataFrame(rows).to_csv(os.path.join(hc.RESULTS, 'dips3_summary.csv'), index=False)


if __name__ == '__main__':
    main([float(a) for a in sys.argv[1:]] or [6.5575, 6.6917, 6.7849])
