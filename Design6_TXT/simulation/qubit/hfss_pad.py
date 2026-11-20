"""Launch-pad reflection: old Design3_2 pad (200 x 200 um) vs the 2026-11-15 pad (400 x 300 um), HFSS driven.

Model (Design6_TXT chip stack: 500 um Si), --b2b (used): BondPad -> 1 mm straight CPW (10 / 6 um) -> mirrored BondPad,
  lumped 50 Ohm ports on both pads' back gaps (pin 'bond', full pad width) = where the bonds land. The pad parameters are
  extracted from this symmetric 2-port (extract()); one launch is then PCB -> bonds -> pad -> matched CPW.
  (Without --b2b: port 2 at the CPW end; its lumped-port mismatch (~ -13 dB at 7 GHz) masks the pad, not used.)
Wire bonds are added afterwards as a series inductance L_b at port 1 (PCB side 50 Ohm):
  S11 seen from the PCB for L_b = 0 ... 1.2 nH, plus the typical cases
  both pads: 2 parallel bonds (L_b ~ 0.55 nH)   [1 mm, 25 um Al, ~100 um pitch].

Usage (qmetal): python driven/hfss_pad.py [--only old|new] [--analyse-only]
Output: results/<date>/pad_s11/{S_old.npy, S_new.npy, f.npy, pad_s11.json, pad_s11.png}
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

SIM = _S('Design6_TXT/simulation')
DESIGN = _S('Design6_TXT'); REPO = _S('.')
sys.path[:0] = [_I('Design6_TXT/layout'), _I('Design5_6_purcell-res/simulation')]

PADS = dict(old=dict(pad_width='200um', pad_height='200um', pad_gap='76.5um', taper_height='122um', back_gap='76.5um'),
            new=None)                                   # new = design6_layout.LP
L_TYP = dict(old=0.55e-9, new=0.55e-9)     # 2 bonds on both pads (user, 2026-11-15)


def build(kind, b2b=False):
    import design6_layout as L
    from bondpad import BondPad
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
    d = L.build_design(channels=False, control_lines=False)
    d.delete_all_components()
    o = dict(L.LP) if kind == 'new' else dict(L.LP, **PADS['old'])
    BondPad(d, 'pad', options=dict(o, pos_x='0mm', pos_y='0mm', orientation='0', **{
        'trace_width': L.CPW_CTL['trace_width'], 'trace_gap': L.CPW_CTL['trace_gap']}))
    if b2b:     # back-to-back: pad -> 1 mm CPW -> mirrored pad, both ports on the bond planes (no CPW-end port)
        BondPad(d, 'pad2', options=dict(o, pos_x='1.06mm', pos_y='0mm', orientation='180', **{
            'trace_width': L.CPW_CTL['trace_width'], 'trace_gap': L.CPW_CTL['trace_gap']}))
        RouteStraight(d, 'cpw', options=dict(pin_inputs=dict(start_pin=dict(component='pad', pin='tie'),
                                                             end_pin=dict(component='pad2', pin='tie')), **L.CPW_CTL))
        return d
    OpenToGround(d, 'end', options=dict(pos_x='1.03mm', pos_y='0mm', orientation='0'))
    RouteStraight(d, 'cpw', options=dict(pin_inputs=dict(start_pin=dict(component='pad', pin='tie'),
                                                         end_pin=dict(component='end', pin='open')), **L.CPW_CTL))
    return d


def run(kind, out, a):
    from qiskit_metal import Dict
    from qiskit_metal.analyses.simulation.scattering_impedance import ScatteringImpedanceSim
    import save_aedt
    save_aedt.use_project(os.environ.get('AEDT_PROJECT', _S('Design6_TXT/simulation/aedt/Design6_flux.aedt')))
    d = build(kind, a.b2b)
    tag = ('b2b_' if a.b2b else '') + kind
    sim = ScatteringImpedanceSim(d, 'hfss')
    for k, v in dict(freq_ghz=8.0, max_delta_s=0.005, max_passes=a.maxp, min_passes=2, min_converged=2,
                     pct_refinement=30, basis_order=1,
                     sweep_setup=Dict(name='Sweep', start_ghz=0.5, stop_ghz=12.0, count=231, step_ghz=None,
                                      type='Interpolating', save_fields=False)).items():
        sim.setup[k] = v
    t0 = time.time()
    comps, ports = ((['pad', 'cpw', 'pad2'], [('pad', 'bond', 50), ('pad2', 'bond', 50)]) if a.b2b else
                    (['pad', 'cpw'], [('pad', 'bond', 50), ('cpw', 'end', 50)]))
    sim.run_sim(name=save_aedt.design_name(f'pad_{tag}'), components=comps,
                open_terminations=[], port_list=ports, box_plus_buffer=True)
    names = ['S11', 'S12', 'S21', 'S22']
    freqs, curves, _ = sim.renderer.get_params(names)
    S = np.zeros((len(freqs), 2, 2), complex)
    for n, c in zip(names, curves):
        S[:, int(n[1]) - 1, int(n[2]) - 1] = np.asarray(c)
    np.save(_relocate(os.path.join(out, f'S_{tag}.npy')), S); np.save(_relocate(os.path.join(out, 'f.npy')), np.asarray(freqs, float))
    print(tag, f'solved in {(time.time() - t0) / 60:.1f} min', flush=True)
    try:
        save_aedt.snapshot(out, f'pad_{tag}')
    except Exception as ex:
        print('snapshot failed:', repr(ex)[:100])


def extract(S, f, z0=50.0):
    """Back-to-back pad - 1 mm CPW - pad (--b2b): pi model of the symmetric 2-port at 2 GHz -> per-pad excess shunt C
    and series L (the 1 mm 10 / 6 um CPW, 50.8 Ohm, sqrt(eps_eff) 2.49, is subtracted), i.e. each pad + taper as a
    short line Z_p = sqrt(L/C), delay sqrt(L C). Port 2 of the HFSS run has reversed polarity: S21, S12 -> -S21, -S12."""
    S = S.copy(); S[:, 0, 1] *= -1; S[:, 1, 0] *= -1
    I = np.eye(2); w = 2 * np.pi * f
    Y = np.array([np.linalg.solve(I + s, I - s) / z0 for s in S])
    i = np.argmin(abs(f - 2e9))
    Csh, Lse = (Y[i, 0, 0] + Y[i, 0, 1]).imag / w[i], (-1 / Y[i, 0, 1]).imag / w[i]
    Lline, Cline = 50.8 * 2.49 / 3e8 * 1e-3, 2.49 / (3e8 * 50.8) * 1e-3
    Cp, Lp = Csh - Cline / 2, (Lse - Lline) / 2
    return Cp, Lp


def launch_s11(Cp, Lp, Lb, f, z0=50.0):
    """PCB (50 Ohm) -> bonds (series Lb) -> pad section (Z_p, delay) -> matched 50 Ohm CPW."""
    Zp, tau = np.sqrt(Lp / Cp), np.sqrt(Lp * Cp)
    t = np.tan(2 * np.pi * f * tau)
    zin = Zp * (z0 + 1j * Zp * t) / (Zp + 1j * z0 * t) + 1j * 2 * np.pi * f * Lb
    return (zin - z0) / (zin + z0)


def analyse(out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f = np.load(_relocate(os.path.join(out, 'f.npy')))
    db = lambda x: 20 * np.log10(np.abs(x))
    res, fig = {}, plt.figure(figsize=(12, 4.4))
    ax = fig.subplots(1, 2)
    lbs = np.linspace(0, 1.2, 61)
    for k in ('old', 'new'):
        fn = _relocate(os.path.join(out, f'S_b2b_{k}.npy'))
        if not os.path.exists(fn):
            continue
        Cp, Lp = extract(np.load(fn), f)
        c = 'C0' if k == 'old' else 'C3'
        res[k] = dict(C_excess_pF=Cp * 1e12, L_pad_nH=Lp * 1e9, Z_pad_ohm=float(np.sqrt(Lp / Cp)),
                      delay_ps=float(np.sqrt(Lp * Cp) * 1e12), L_bond_compensating_nH=Cp * 50 ** 2 * 1e9,
                      S11_dB={f'Lb {lb} nH': {f'{g} GHz': round(float(db(launch_s11(Cp, Lp, lb * 1e-9, g * 1e9))), 1)
                                              for g in (4, 5, 6, 7, 8)} for lb in (0, 0.3, 0.55, 0.8, 1.0)},
                      typical_L_bond_nH=L_TYP[k] * 1e9)
        ax[0].plot(f / 1e9, db(launch_s11(Cp, Lp, L_TYP[k], f)), c, lw=1.5,
                   label=f'{k} pad, 2 bonds ({L_TYP[k] * 1e9:.2f} nH)')
        ax[0].plot(f / 1e9, db(launch_s11(Cp, Lp, 1e-9, f)), c + ':', lw=1.2, label=f'{k} pad, 1 bond (1 nH)')
        ax[1].plot(lbs, [float(db(launch_s11(Cp, Lp, lb * 1e-9, 7e9))) for lb in lbs], c, label=f'{k} pad')
        ax[1].axvline(L_TYP[k] * 1e9, c=c, ls=':', lw=1)
    ax[0].axvspan(4, 7.2, color='0.9', zorder=0)
    ax[0].set_xlabel('f (GHz)'); ax[0].set_ylabel('|S11| (dB)'); ax[0].set_ylim(-45, 0); ax[0].legend(fontsize=8)
    ax[0].set_title('one launch, PCB 50 Ohm -> bonds -> pad -> CPW')
    ax[1].set_xlabel('bond inductance L_b (nH)'); ax[1].set_ylabel('|S11| at 7 GHz (dB)'); ax[1].legend()
    ax[1].set_title('dotted: 2 bonds (0.55 nH)')
    fig.tight_layout(); fig.savefig(_relocate(os.path.join(out, 'pad_s11.png')), dpi=130); plt.close(fig)
    json.dump(res, open(_relocate(os.path.join(out, 'pad_s11.json')), 'w'), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', choices=['old', 'new'])
    ap.add_argument('--maxp', type=int, default=10)
    ap.add_argument('--analyse-only', action='store_true')
    ap.add_argument('--b2b', action='store_true', help='back-to-back pad - CPW - pad, ports on both bond planes')
    ap.add_argument('--dry', action='store_true', help='build the two test designs only (no HFSS)')
    ap.add_argument('--date', default=time.strftime('%Y-%m-%d'))
    a = ap.parse_args()
    out = _relocate(os.path.join(SIM, 'results', a.date, 'pad_s11')); os.makedirs(out, exist_ok=True)
    if a.dry:
        for k in ('old', 'new'):
            d = build(k)
            print(k, {n: (np.round(c.pins[p].middle, 4).tolist(), round(c.pins[p].width * 1e3, 1), round(c.pins[p].gap * 1e3, 1))
                      for n, c in d.components.items() for p in c.pins})
        sys.exit()
    if not a.analyse_only:
        for k in ([a.only] if a.only else ['old', 'new']):
            run(k, out, a)
    analyse(out)
