"""Charging energies and couplings of one FTF set from a Maxwell capacitance matrix.

Method (see capacitance_method.md, thesis eqs. A.14 - A.15):

    nodes  = [FX1 pad a, FX1 pad b, coupler, FX2 pad a, FX2 pad b, (readout 1, readout 2)]
    C~     = M^-T C M^-1      (M: sum / difference coordinates of each FX pad pair)
    Ec_i   = e^2/2 * C~^-1[i, i]                    i = 1 (FX1), 2 (coupler), 4 (FX2)
    J_ij   = 4 e^2 * C~^-1[i, j]                    (1,2) (2,4) (1,4)
    readout g: readout_g() (beta * V_zpf, lambda/2 resonator)
    charge line: charge_line()

Every conductor that is not in `nodes` (ground, flux lines, charge lines, ...) is
treated as grounded: its rows / columns are simply dropped from the Maxwell matrix.
Indices above are 0-based, as in the thesis.
"""
import numpy as np
import pandas as pd
from scipy.constants import e, h, hbar

# Q3D net names produced by qiskit-metal: '<polygon name>_<component name>'
DEFAULT_NODES = [
    'pad_bot_FX1', 'pad_top_FX1',      # 0, 1  -> FX1 sum / difference
    'tmon_TC',                         # 2     -> coupler
    'pad_top_FX2', 'pad_bot_FX2',      # 3, 4  -> FX2 sum / difference
    'readout_line_FX1',                # 5     -> resonator 1
    'readout_line_FX2',                # 6     -> resonator 2
]


def transform_matrix(n, pairs=((0, 1), (3, 4))):
    """M of eq. A.14: rows (1, 1) and (1, -1) for each pad pair, identity elsewhere."""
    M = np.eye(n)
    for a, b in pairs:
        M[a, a], M[a, b] = 1, 1
        M[b, a], M[b, b] = 1, -1
    return M


def reduced_matrix(cmat, nodes=DEFAULT_NODES):
    """Maxwell matrix (fF, DataFrame from Q3D) restricted to `nodes`, in farad."""
    missing = [n for n in nodes if n not in cmat.index]
    if missing:
        raise KeyError(f'nodes {missing} not in the capacitance matrix: {list(cmat.index)}')
    return cmat.reindex(index=nodes, columns=nodes).to_numpy(dtype=float) * 1e-15


def analyse(cmat, nodes=DEFAULT_NODES, verbose=True):
    """Return a dict with Ec (GHz) and J (MHz) of the FTF set."""
    C = reduced_matrix(cmat, nodes)
    M = transform_matrix(len(nodes))
    Mi = np.linalg.inv(M)
    Ct = Mi.T @ C @ Mi                       # C~ = M^-T C M^-1
    Ci = np.linalg.inv(Ct)

    Ec = lambda i: e**2 / 2 * Ci[i, i] / h / 1e9
    J = lambda i, j: 4 * e**2 * Ci[i, j] / h / 1e6

    res = {
        'Ec_FX1 [GHz]': Ec(1),
        'Ec_TC [GHz]': Ec(2),
        'Ec_FX2 [GHz]': Ec(4),
        'J_1c [MHz]': J(1, 2),
        'J_2c [MHz]': J(2, 4),
        'J_12 [MHz]': J(1, 4),
    }
    if verbose:
        print(pd.Series(res).round(4).to_string())
    return res


def total_capacitances(cmat, nodes=DEFAULT_NODES):
    """C_sigma (fF) seen by each mode = 1 / C~^-1[i, i]; handy for tuning."""
    C = reduced_matrix(cmat, nodes)
    Mi = np.linalg.inv(transform_matrix(len(nodes)))
    Ci = np.linalg.inv(Mi.T @ C @ Mi)
    return {'FX1': 1e15 / Ci[1, 1], 'TC': 1e15 / Ci[2, 2], 'FX2': 1e15 / Ci[4, 4]}


def target_C(Ec_GHz):
    """C_sigma (fF) that gives a charging energy Ec (GHz)."""
    return e**2 / (2 * h * Ec_GHz * 1e9) * 1e15


def find_node(cmat, *keys):
    """First net name containing any of `keys` (Q3D names a merged net after one of its objects)."""
    for k in keys:
        for n in cmat.index:
            if k in n:
                return n
    raise KeyError(f'no net matching {keys} in {list(cmat.index)}')


def readout_g(cmat, fx_nodes, res_node, Z0=45.0, f_r=None, lam=2, transmon_EJ_over_Ec=None, verbose=True):
    """Readout coupling of one fluxonium to its resonator (resonator rendered in Q3D).

    fx_nodes = (pad a, pad b) of a fluxonium, or (island,) of a grounded transmon;
    res_node = net of resonator + readout pad.
    All other conductors are treated as grounded (coupler, other FX: small effect on g).

        beta  = C~^-1[q, r] / C~^-1[r, r]          voltage-division ratio (= C_c / C_sigma,q)
        V_zpf = sqrt(hbar w_r / 2 C_eff)          lambda/2: C_eff = pi / (2 w_r Z0) -> V_zpf = w_r sqrt(hbar Z0 / pi)
                                                  lambda/4: C_eff = pi / (4 w_r Z0)
        hbar g = 2 e beta V_zpf                   (thesis convention: per charge number n,
                                                   fluxonium matrix element excluded)
    transmon_EJ_over_Ec: also return the transmon-convention g (A.21), g * (32 Ec/EJ)^(1/4).
    Pass f_r (Hz, the resonator design frequency). Without it, f_r is estimated as
    1 / (2 Z0 C_line) from the Q3D self capacitance of the resonator net, which comes
    out too high when the tee section is not rendered / the mesh is coarse.
    """
    nodes = list(fx_nodes) + [res_node]
    C = reduced_matrix(cmat, nodes)
    if len(fx_nodes) == 2:                       # floating qubit: difference coordinate
        Mi = np.linalg.inv(transform_matrix(3, pairs=((0, 1),)))
        Ci = np.linalg.inv(Mi.T @ C @ Mi)
        q, r = 1, 2
    else:                                        # grounded qubit (transmon coupler): one node
        Ci = np.linalg.inv(C)
        q, r = 0, 1
    beta = Ci[q, r] / Ci[r, r]
    C_line = cmat.loc[res_node, res_node] * 1e-15            # total capacitance of the resonator net
    if f_r is None:
        f_r = 1 / (2 * Z0 * C_line) if lam == 2 else 1 / (4 * Z0 * C_line)
    w_r = 2 * np.pi * f_r
    C_eff = np.pi / (2 * w_r * Z0) if lam == 2 else np.pi / (4 * w_r * Z0)
    V_zpf = np.sqrt(hbar * w_r / (2 * C_eff))
    g = 2 * e * beta * V_zpf / hbar / (2 * np.pi)
    out = {'beta = Cc/Csigma': beta, 'C_line [fF]': C_line * 1e15, 'f_r [GHz]': f_r / 1e9,
           'V_zpf [uV]': V_zpf * 1e6, 'g/2pi [MHz]': g / 1e6}
    if transmon_EJ_over_Ec:
        # thesis A.21: transmon convention multiplies by (32 Ec / EJ)^(1/4)
        out['g_transmon/2pi [MHz]'] = g / 1e6 * (32 / transmon_EJ_over_Ec) ** 0.25
    if verbose:
        print(f'  {res_node}:', ', '.join(f'{k} = {v:.4g}' for k, v in out.items()))
    return out


def charge_line(cmat, qubit_nodes, line_node, f_q, n_q, R=50.0, rabi=10e6):
    """Charge-line drive of a qubit through a 50 Ohm line, from the Maxwell matrix.

    beta = C~^-1[q, cl] / C~^-1[cl, cl] (as for readout). With a charge matrix element n_q
    at frequency f_q:
        T1 limit      Gamma = 8 e^2 beta^2 |n_q|^2 w R / hbar     (T = 0 noise of R)
        drive         V for a Rabi rate `rabi`:  V = h * rabi / (2 e beta |n_q|)
    """
    nodes = list(qubit_nodes) + [line_node]
    C = reduced_matrix(cmat, nodes)
    if len(qubit_nodes) == 2:
        Mi = np.linalg.inv(transform_matrix(3, pairs=((0, 1),)))
        Ci = np.linalg.inv(Mi.T @ C @ Mi)
        q, r = 1, 2
    else:
        Ci = np.linalg.inv(C)
        q, r = 0, 1
    beta = Ci[q, r] / Ci[r, r]
    w = 2 * np.pi * f_q
    gamma = 8 * e**2 * beta**2 * n_q**2 * w * R / hbar
    V = h * rabi / (2 * e * abs(beta) * n_q)
    P_dbm = 10 * np.log10(V**2 / (2 * R) / 1e-3)
    C_direct = -sum(cmat.loc[line_node, n] for n in qubit_nodes)
    return {'C_line-qubit [fF]': C_direct, 'beta': beta, 'T1 limit [us]': 1e6 / gamma,
            f'P for {rabi / 1e6:.0f} MHz Rabi [dBm at chip]': P_dbm}
