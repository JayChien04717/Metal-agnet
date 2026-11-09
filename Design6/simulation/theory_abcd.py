"""Circuit (ABCD / admittance) model of the Design6 channel using the Q3D capacitances, and an HFSS-vs-theory
comparison per run.

    P1 - TL - [C_in pi] - TL - node n0 -- TL - P2
                           |  n0 shunt branch: TL - [C_p pi] - filter TL (open at C_p, shorted at the far end)
                                                              | tap at D_CJ from C_p: [C_j pi] - resonator TL (open end at d_r, short at L_res)
Each enclosing capacitor / C_j is a pi network: series C (Q3D coupling) + the two electrode capacitances to ground
(Q3D: C_ii + C_ij). Lines: Z0 and v of the CPW 10/6 on Si, lossless.

Usage:  python theory_abcd.py <tag>      -> results/<run folder>/<tag>_theory_vs_hfss.png + theory csv
        archive(tag) is called at the end of hfss_chan6.py
"""
import os, sys, json, time, shutil, glob, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, 'results')
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'pattern design')); sys.path.insert(1, HERE)


def q3d_caps(name):
    """(series C, ground C of electrode 1, ground C of electrode 2) in fF from the latest Q3D matrix."""
    cm = pd.read_csv(os.path.join(RES, f'q3d6_{name}_cmat_fF.csv'), index_col=0)
    sig = [c for c in cm.index if 'ground' not in c.lower()]
    ser = -cm.loc[sig[0], sig[1]]
    return ser, cm.loc[sig[0], sig[0]] - ser, cm.loc[sig[1], sig[1]] - ser, sig


def load_block(path_csv, zref=50.0):
    """HFSS 3-port of the C_j branch (hfss_block6.py) -> dict(f, Y) with Y the (n, 3, 3) admittance matrix."""
    S = pd.read_csv(path_csv, index_col=0).map(complex)
    f = S.index.values.astype(float)
    Sm = np.stack([[S[f'S{i}{j}'].values for j in (1, 2, 3)] for i in (1, 2, 3)]).transpose(2, 0, 1)
    I = np.eye(3)
    Z = zref * np.linalg.solve(I - Sm, I + Sm)            # Z = z0 (I - S)^-1 (I + S)
    return dict(f=f, Y=np.linalg.inv(Z), path=path_csv)


def block_Yin(f_ghz, blk, YL2, YL3):
    """Input admittance at port 1 (C_p pin b) with ports 2 / 3 loaded by YL2 / YL3:
    Y_in = Y11 - Y1x (Yxx + diag(YL))^-1 Yx1  (node equations, load current I_k = -YL_k V_k)."""
    f = np.real(np.asarray(f_ghz, complex)).ravel()
    Yb = np.empty((len(f), 3, 3), complex)
    for i in range(3):
        for j in range(3):
            Yb[:, i, j] = np.interp(f, blk['f'], blk['Y'][:, i, j].real) + 1j * np.interp(f, blk['f'], blk['Y'][:, i, j].imag)
    YL2 = np.broadcast_to(YL2, f.shape); YL3 = np.broadcast_to(YL3, f.shape)
    Yxx = Yb[:, 1:, 1:].copy(); Yxx[:, 0, 0] += YL2; Yxx[:, 1, 1] += YL3
    Vx = np.linalg.solve(Yxx, Yb[:, 1:, 0][..., None])[..., 0]
    out = Yb[:, 0, 0] - np.einsum('ni,ni->n', Yb[:, 0, 1:], Vx)
    return out.reshape(np.shape(f_ghz))


def branch_Y(f_ghz, P):
    """Admittance of the shunt branch (stub - C_p - filter - C_j - resonator) seen at node n0."""
    Z0, v = P['Z0'], P['v']
    w = 2 * np.pi * np.asarray(f_ghz) * 1e9
    beta = w / v
    Y0 = 1 / Z0
    t = lambda l: np.tan(beta * l * 1e-3)
    Yopen = lambda l: 1j * Y0 * t(l)
    Yshort = lambda l: -1j * Y0 / t(l)
    cjs, cjg1, cjg2 = P['cj']
    Yres = 1j * w * cjg2 * 1e-15 + Yopen(P['d_r']) + Yshort(P['L_res'] - P['d_r'])      # resonator node
    Yser = (1j * w * cjs * 1e-15 * Yres) / (1j * w * cjs * 1e-15 + Yres)
    Ytap = 1j * w * cjg1 * 1e-15 + Yser + Yshort(P['L_filt'] - P['D_cj'])                # filter node at the tap
    # transfer back to the C_p pin b along the open filter line of length D_cj
    td = t(P['D_cj'])
    Yf = Y0 * (Ytap + 1j * Y0 * td) / (Y0 + 1j * Ytap * td)
    if P.get('block') is not None:                       # hybrid: HFSS 3-port of the C_j branch replaces the ideal one
        Yf = block_Yin(f_ghz, P['block'], Yshort(P['Lf_rest']), Yshort(P['Lr_rest']))
    # C_p pi: frame side (pin a, toward n0) / plate side (pin b, toward filter)
    cps, cpg_a, cpg_b = P['cp']
    Ycp = Yf
    Y1 = 1j * w * cpg_b * 1e-15 + Ycp
    Y2 = (1j * w * cps * 1e-15 * Y1) / (1j * w * cps * 1e-15 + Y1)
    Y3 = 1j * w * cpg_a * 1e-15 + Y2
    # stub line (n0 -> C_p pin a): input admittance of a line of length l_stub loaded with Y3
    ts = t(P['l_stub'])
    return Y0 * (Y3 + 1j * Y0 * ts) / (Y0 + 1j * Y3 * ts)


def s21_theory(f_ghz, P):
    """P: dict with Z0, v (m/s), caps and lengths (mm). Returns complex S11, S21."""
    Z0, v = P['Z0'], P['v']
    w = 2 * np.pi * np.asarray(f_ghz) * 1e9
    beta = w / v
    Y0 = 1 / Z0
    tl = lambda l: np.array([[np.cos(beta * l * 1e-3), 1j * Z0 * np.sin(beta * l * 1e-3)],
                             [1j * Y0 * np.sin(beta * l * 1e-3), np.cos(beta * l * 1e-3)]])
    # vectorised 2x2 over frequency: arrays of shape (2, 2, n)
    def mul(A, B):
        return np.einsum('ijn,jkn->ikn', A, B)
    one = np.zeros((2, 2, len(w)), complex); one[0, 0] = 1; one[1, 1] = 1

    def shunt(Y):
        M = one.copy(); M[1, 0] = Y; return M

    def series(Z):
        M = one.copy(); M[0, 1] = Z; return M

    def pi(cser, cga, cgb):
        return mul(mul(shunt(1j * w * cga * 1e-15), series(1 / (1j * w * cser * 1e-15))), shunt(1j * w * cgb * 1e-15))

    Ybr = branch_Y(f_ghz, P)
    cis, cig_a, cig_b = P['cin']
    M = mul(tl(P['l_feed_in']), pi(cis, cig_a, cig_b))
    M = mul(M, tl(P['l_feed_c']))
    M = mul(M, shunt(Ybr))
    M = mul(M, tl(P['l_out']))
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    den = A + B / Z0 + C * Z0 + D
    return (A + B / Z0 - C * Z0 - D) / den, 2 / den


def params_from_layout(meta=None):
    import design6_layout as d6, ftf_layout as f3
    d = d6.build_design()
    P = dict(Z0=d6.Z0, v=d6.V)
    for n in ('cin', 'cp', 'cj'):
        s, g1, g2, sig = q3d_caps(n)
        P[n] = (s, g1, g2)
    # EnclosingCap electrodes: the 'frame' net is on pin a, the plate ('lead_b') on pin b (T-type: a = line side).
    # Q3D net order is frame, ..., lead_b; q3d_caps returns [frame, lead_b]
    P['cin'] = P['cin'] if d6.CIN_TYPE == 'T' else (P['cin'][0], P['cin'][2], P['cin'][1])
    P['cp'] = P['cp'] if d6.CP_TYPE == 'T' else (P['cp'][0], P['cp'][2], P['cp'][1])
    pl = lambda n: f3.path_length(d, n)
    P['l_feed_in'] = pl('feed_in') + 0.04            # + lead of the cap
    P['l_feed_c'] = pl('feed_c') + 0.08 + 0.15       # leads + half of the n0 prime line
    P['l_out'] = pl('out_line') + 0.15
    P['l_stub'] = 0.1 + pl('stub_p') + 0.04
    P['L_filt'] = meta['L_FILTER'] if meta else d6.L_FILTER
    P['D_cj'] = d6.D_CJ
    P['L_res'] = meta['L_RES_TOTAL'] if meta else d6.L_RES_TOTAL
    P['d_r'] = d6.VERT + d6.OPEN_DY - d6.D_CJ
    return P


def archive(tag, run_dir=None):
    """Move this run's files into results/<YYYYmmdd_HHMM>_<tag>/ and add the theory curve and the overlay plot."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    stamp = time.strftime('%Y%m%d_%H%M')
    run_dir = run_dir or os.path.join(RES, f'{stamp}_{tag}')
    os.makedirs(run_dir, exist_ok=True)
    for p in glob.glob(os.path.join(RES, f'{tag}_*')):
        shutil.move(p, os.path.join(run_dir, os.path.basename(p)))
    for n in ('cin', 'cp', 'cj'):                      # the capacitances used by the theory curve
        src = os.path.join(RES, f'q3d6_{n}_cmat_fF.csv')
        if os.path.exists(src):
            shutil.copy(src, os.path.join(run_dir, f'{tag}_q3d_{n}_cmat_fF.csv'))
    db = pd.read_csv(os.path.join(run_dir, f'{tag}_dB.csv'), index_col=0)
    P = params_from_layout(json.load(open(os.path.join(run_dir, f'{tag}_meta.json'))))
    f = db.index.values
    s11, s21 = s21_theory(f, P)
    th = pd.DataFrame({'S11_dB': 20 * np.log10(np.abs(s11)), 'S21_dB': 20 * np.log10(np.abs(s21))}, index=db.index)
    th.to_csv(os.path.join(run_dir, f'{tag}_theory_dB.csv'))
    json.dump({k: (list(map(float, v)) if isinstance(v, tuple) else v) for k, v in P.items()},
              open(os.path.join(run_dir, f'{tag}_theory_params.json'), 'w'), indent=1)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    a1.plot(f, db.S21_dB, label='HFSS S21'); a1.plot(f, th.S21_dB.clip(lower=-90), '--', label='theory S21 (ABCD, Q3D C)')
    a2.plot(f, db.S11_dB, label='HFSS S11'); a2.plot(f, th.S11_dB, '--', label='theory S11')
    a1.set_ylabel('S21 (dB)'); a2.set_ylabel('S11 (dB)'); a2.set_xlabel('f (GHz)')
    for a in (a1, a2):
        a.grid(alpha=.3); a.legend()
    a1.set_title(f'{tag}  ({stamp})   Cin {P["cin"][0]:.1f}, Cp {P["cp"][0]:.1f}, Cj {P["cj"][0]:.2f} fF')
    fig.tight_layout(); fig.savefig(os.path.join(run_dir, f'{tag}_theory_vs_hfss.png'), dpi=120); plt.close(fig)
    return run_dir


if __name__ == '__main__':
    print(archive(sys.argv[1]))
