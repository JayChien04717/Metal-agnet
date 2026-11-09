"""ABCD check of the kappa=3 MHz design: plain line (3 hangers) and with Purcell caps, vs HFSS."""
import os, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import abcd_model as am
import design3_purcell as dp, purcell_layout as pl
from scipy.signal import find_peaks
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

names = list(dp.RESONATORS)
d = dp.build_design()
# node positions as in abcd_model (centre of each tee prime line along filter line)
xs = {n: float(str(d.components[f'tee_{n}'].options.pos_x).replace('mm', '')) for n in names}
x_first = min(xs.values())
gaps = [xs[b] - xs[a] for a, b in zip(names[:-1], names[1:])]
res = dp.RESONATORS
def hang(n): return am.hanger(res[n]['Cc'], res[n]['length'])
# plain line: hangers separated by the tee-to-tee distance (abs x spacing, mirror ignored)
pos = [abs(xs[names[0]]), ]
plain = [hang(names[0])]
for g, n in zip(gaps, names[1:]):
    plain += [am.tl(abs(g)), hang(n)]
Sp = am.s21(plain); sp = am.db(Sp)
idx, pr = find_peaks(-sp, prominence=1.0)
print('ABCD plain line dips:', [(round(am.F[i]/1e9, 4), round(sp[i], 1)) for i in idx])
print('HFSS plain dips: 6.585, 6.749, 6.800, 6.899')
print('design f:', {n: res[n]['f'] for n in names}, 'length', {n: res[n]['length'] for n in names})

# with caps: filter path PIN_LENGTH, hangers at tee centres
cap = pl.CAP_FF
delta = 0.98
tot = pl.PIN_LENGTH
mlen = pl.path_length(d, 'filter_in')
x1 = xs[names[0]] - dp.RESONATORS[names[0]]['coupling_length']
nodes = [mlen + (xs[n] - x1) for n in names]
chain = [am.series_c(cap), am.tl(delta/2 + nodes[0])]
for i, n in enumerate(names):
    chain.append(hang(n))
    chain.append(am.tl((nodes[i+1]-nodes[i]) if i < len(names)-1 else tot - nodes[-1] + delta/2))
chain.append(am.series_c(cap))
Sc = am.s21(chain); sc = am.db(Sc)
filt = am.db(am.s21([am.series_c(cap), am.tl(tot + delta), am.series_c(cap)]))
print('filter alone peak (GHz):', am.F[np.argmax(filt)]/1e9)
idx, pr = find_peaks(-sc, prominence=1.0)
print('ABCD with Purcell dips:', [(round(am.F[i]/1e9, 4), round(sc[i], 1)) for i in idx])
print('ABCD passband 3dB edges of filter alone:', am.F[filt > filt.max()-3][[0, -1]]/1e9)
fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
ax[0].plot(am.F/1e9, sp); ax[0].set_title('ABCD plain line'); 
ax[1].plot(am.F/1e9, sc, label='with Purcell'); ax[1].plot(am.F/1e9, filt, '--', label='filter alone'); ax[1].legend()
for a in ax:
    a.set_ylim(-60, 2); a.grid(alpha=.3)
    for n in names: a.axvline(res[n]['f'], color='gray', ls=':')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'results', 'abcd_check.png'), dpi=110)
