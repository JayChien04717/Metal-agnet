"""Plot both channels of a loaded HFSS sweep with its complex four-pole fit."""
from argparse import ArgumentParser
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from compare_loaded_trim_r9 import join_grids
from compare_loaded_s21_r9 import db, read
from readout_poles import load


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('solver_run', type=Path)
    parser.add_argument('fit_run', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if str(args.solver_run.resolve()) not in read(args.fit_run / 'run.json')['upstream']:
        parser.error('The fit must use the selected HFSS run')

    wide = load(args.solver_run / 'raw/S_wide.csv')
    with np.load(args.fit_run / 'raw/fitted_s21.npz') as fit:
        f_fit, s_fit = fit['f_GHz'], fit['S21_fit']
        freq, s21 = join_grids(wide, (f_fit, fit['S21']))
        poles = fit['poles_GHz'].real

    fig, ax = plt.subplots(figsize=(10.8, 4.5), layout='constrained')
    ax.plot(freq, db(s21), color='#193e67', lw=1.8, label='HFSS $S_{21}$')
    ax.plot(f_fit, db(s_fit), color='#df642e', lw=1.5, ls='--', label='Four-pole fit')
    for label, pair in (('Q2', poles[:2]), ('Q1', poles[2:])):
        ax.text(float(pair.mean()), 1.4, label, ha='center', va='center', fontweight='semibold')
    ax.set(xlim=(freq[0], freq[-1]), ylim=(min(db(s21)) - 3, 3.7),
           xlabel='Frequency (GHz)', ylabel='$|S_{21}|$ (dB)')
    ax.grid(alpha=.22)
    ax.legend(loc='lower left')
    fig.savefig(args.output, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
